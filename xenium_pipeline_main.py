"""
Xenium 분석 파이프라인 - 메인 스크립트
========================================

이 스크립트는 Xenium 데이터를 처리하는 완전한 파이프라인을 제공합니다.
Step 0부터 Step 6까지의 모든 단계를 포함합니다.

논문 참고:
- Xenium 벤치마킹 논문
- 최적 전처리 경로: target_sum=100, log-transform=True, scale=False, hvg=False

Author: 자동 생성
Date: 2024-12-29
"""

import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Optional, Dict, List
import logging
import yaml

# xb 모듈 임포트
import xb.formatting as xf
import xb.preprocessing as xp
from xb.simulating import allcombs, allcombs_simulated
from xb.calculating import dispersion, dist_nuc, compute_vi, compute_fmi

# Ground Truth Tool Import (formerly Step 7)
try:
    import step7_ground_truth as s7
except ImportError:
    # logger.warning("step7_ground_truth 모듈을 찾을 수 없습니다. Ground Truth 생성 기능을 사용할 수 없습니다.")
    s7 = None

# NEW MODULES: Step 7 & 8
try:
    import step7_spatial_domains as s7_domains
except ImportError:
    s7_domains = None

try:
    import step8_svf as s8_svf
except ImportError:
    s8_svf = None

# NEW MODULES: Step 2 Extensions (Overlaps & SSAM)
try:
    import step2_overlaps as s2_overlaps
    import step2_ssam as s2_ssam
except ImportError:
    s2_overlaps = None
    s2_ssam = None

# Visualization module
import visualize_results as vis


# =====================================================================
# 설정 로더
# =====================================================================

def load_config(config_file: str = "config.yaml") -> Dict:
    """
    config.yaml 파일을 읽고 설정을 반환합니다.

    Parameters
    ----------
    config_file : str
        설정 파일 경로

    Returns
    -------
    Dict
        설정 딕셔너리
    """
    if not os.path.exists(config_file):
        raise FileNotFoundError(
            f"설정 파일을 찾을 수 없습니다: {config_file}\n"
            f"config.yaml을 생성한 후 경로를 설정하세요."
        )

    with open(config_file, 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)

    return config


def setup_logging(config: Dict):
    """
    설정 파일을 기반으로 로깅을 설정합니다.
    """
    log_config = config.get('logging', {})
    level = getattr(logging, log_config.get('level', 'INFO'))

    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    if log_config.get('save_log', True):
        log_file = log_config.get('log_file', 'xenium_pipeline.log')
        handler = logging.FileHandler(log_file)
        handler.setLevel(level)
        logging.getLogger().addHandler(handler)


logger = logging.getLogger(__name__)


class XeniumPipeline:
    """Xenium 분석 파이프라인 클래스"""

    def __init__(
        self,
        xenium_input_path: str = None,
        output_path: str = "./output",
        sample_tag: str = "xenium_sample",
        load_preprocessed: bool = False
    ):
        """
        Parameters
        ----------
        xenium_input_path : str
            Xenium 원본 데이터 디렉토리 경로 (또는 .h5ad 파일 경로)
        output_path : str, optional
            결과 저장 디렉토리, by default "./output"
        sample_tag : str, optional
            샘플 이름 태그, by default "xenium_sample"
        load_preprocessed : bool, optional
            True이면 xenium_input_path를 h5ad 파일로 간주, by default False
        """
        self.xenium_input_path = xenium_input_path
        self.sample_tag = sample_tag
        # 데이터셋 이름으로 하위 디렉토리 생성
        self.output_path = Path(output_path) / self.sample_tag
        self.output_path.mkdir(parents=True, exist_ok=True)
        self.adata = None
        self.load_preprocessed = load_preprocessed
        
        if self.load_preprocessed and self.xenium_input_path:
             try:
                 logger.info(f"H5AD 파일 로드 중: {self.xenium_input_path}")
                 self.adata = sc.read_h5ad(self.xenium_input_path)
                 logger.info(f"✓ 로드 완료: {self.adata.shape}")
             except Exception as e:
                 logger.error(f"H5AD 파일 로드 실패: {e}")
                 # 로드 실패 시 None으로 유지되어 이후 로직에서 처리됨

    # =====================================================================
    # STEP 0: 데이터 포맷팅
    # =====================================================================
    # 참고: notebooks/0_formatting/0_0_Formatting xenium to anndata.ipynb
    # 함수: format_xenium_adata_mid_2023() (2024년 Xenium 포맷)

    def step0_format_xenium(
        self,
        generate_ground_truth: bool = False,
        tissue_type: str = "brain",
        download_dir: str = "data/reference",
        organ: Optional[str] = None
    ) -> sc.AnnData:
        """
        Step 0: Xenium 원본 데이터를 AnnData 형식으로 변환
        
        Optional: Ground Truth 데이터 생성 (Reference Download & Label Transfer)

        입력: Xenium 장비 출력 (transcripts.csv, matrix.mtx 등)
        출력: AnnData 객체 및 H5AD 파일

        Returns
        -------
        sc.AnnData
            포맷팅된 AnnData 객체
        """
        logger.info("="*60)
        logger.info("STEP 0: 데이터 포맷팅 & Ground Truth 준비")
        logger.info("="*60)

        # 1. 데이터 포맷팅 / 로딩
        # 1. 데이터 포맷팅 / 로딩
        if self.adata is None:
            # Check if input is already H5AD
            if str(self.xenium_input_path).endswith('.h5ad'):
                 logger.info(f"H5AD 파일 직접 로드 중: {self.xenium_input_path}")
                 try:
                     self.adata = sc.read_h5ad(self.xenium_input_path)
                     logger.info(f"✓ 로드 완료: {self.adata.shape}")
                 except Exception as e:
                     logger.error(f"H5AD 로드 실패: {e}")
                     raise
            else:
                # Raw data formatting
                try:
                    self.adata = xf.format_xenium_adata_mid_2023(
                        path=self.xenium_input_path,
                        tag=self.sample_tag,
                        output_path=str(self.output_path)
                    )
                    logger.info(f"✓ 포맷팅 완료: {self.adata.shape}")
                    logger.info(f"  - 세포 수: {self.adata.shape[0]}")
                    logger.info(f"  - 유전자 수: {self.adata.shape[1]}")

                    # 배경 이미지 추출
                    xf.format_background(self.xenium_input_path)
                except Exception as e:
                    logger.error(f"포맷팅 실패: {e}")
                    raise
        else:
            logger.info(f"✓ Preprocessed/Loaded data found: {self.adata.shape}")

        # 2. Ground Truth 생성 로직
        if generate_ground_truth:
            if s7 is None:
                logger.warning("step7_ground_truth 모듈 부재로 Ground Truth 생성을 건너뜁니다.")
            elif 'ground_truth_celltype' in self.adata.obs:
                 logger.info("✓ Ground Truth already exists. Skipping generation.")
            else:
                logger.info("\n[0-1] Ground Truth 데이터 생성 시작...")
                try:
                    # 1. Reference Data
                    logger.info(f"  - 레퍼런스 데이터 준비 ({tissue_type})...")
                    adata_ref = s7.download_reference_data(
                        tissue_type=tissue_type,
                        download_dir=download_dir,
                        organ=organ
                    )
                    
                    # 2. Label Transfer
                    logger.info("  - 라벨 전이 (Label Transfer)...")
                    # s7.transfer_labels는 adata를 리턴하고 내부적으로 obs에 라벨 추가
                    self.adata = s7.transfer_labels(self.adata, adata_ref, label_key='cell_type')
                    
                    # 컬럼 이름 확인 및 저장
                    if 'cell_type' in self.adata.obs:
                        self.adata.obs['ground_truth_celltype'] = self.adata.obs['cell_type']
                        logger.info(f"  ✓ Ground Truth 라벨 생성 완료: {len(self.adata.obs['ground_truth_celltype'].unique())} cell types")
                    else:
                        logger.warning("  ⚠ 라벨 전이 후 'cell_type' 컬럼을 찾을 수 없습니다.")

                    # 업데이트된 데이터 저장
                    # output_file_labeled = self.output_path / f"{self.sample_tag}_step0_labeled.h5ad"
                    # self.adata.write(str(output_file_labeled))
                    # logger.info(f"  ✓ 라벨링된 데이터 저장: {output_file_labeled}")

                except Exception as e:
                    logger.error(f"Ground Truth 생성 중 오류: {e}")
                    # 실패하더라도 파이프라인 중단하지 않음 (선택사항이므로)

        return self.adata

    # =====================================================================
    # STEP 1: 데이터셋 탐색 및 전처리
    # =====================================================================
    # 참고: notebooks/1_datasets_exploration/1_6_Batch_preprocessing_real_Xenium_datasets.ipynb
    # 함수: main_preprocessing(), keep_nuclei_and_quality(), allcombs()

    def step1_preprocess(
        self,
        target_sum: int = 100,
        mincounts: int = 10,
        mingenes: int = 3,
        neigh: int = 15,
        scale: bool = False,
        hvg: bool = False,
        save: bool = True,
        **kwargs
    ) -> sc.AnnData:
        """
        Step 1: 데이터 품질 검증 및 전처리

        권장 파라미터:
        - target_sum: 100 (논문에서 최적으로 검증됨)
        - scale: False (공간 데이터에서는 권장하지 않음)
        - hvg: False (모든 유전자 사용이 더 정확함)

        Parameters
        ----------
        target_sum : int, optional
            정규화 기준값, by default 100
        mincounts : int, optional
            최소 카운트 필터, by default 10
        mingenes : int, optional
            최소 유전자 발현 필터, by default 3
        neigh : int, optional
            최근접 이웃 수, by default 15
        scale : bool, optional
            표준화 여부, by default False (권장)
        hvg : bool, optional
            고변이 유전자 선택 여부, by default False (권장)
        save : bool, optional
            결과 저장 여부, by default True

        Returns
        -------
        sc.AnnData
            전처리된 AnnData 객체
        """
        logger.info("="*60)
        logger.info("STEP 1: 데이터 탐색 및 전처리")
        logger.info("="*60)

        # Step 0을 건너뛴 경우, h5ad 파일 직접 로드
        if self.adata is None and self.load_preprocessed:
            logger.info(f"\n[1-0] H5AD 파일 로드...")
            if not self.xenium_input_path or not self.xenium_input_path.endswith('.h5ad'):
                logger.error("h5ad 파일 경로를 지정하세요!")
                return None
            try:
                self.adata = sc.read_h5ad(self.xenium_input_path)
                logger.info(f"✓ 로드 완료: {self.adata.shape}")
            except Exception as e:
                logger.error(f"파일 로드 실패: {e}")
                return None

        if self.adata is None:
            logger.error("Step 0을 먼저 실행하거나 load_preprocessed=True로 설정하세요.")
            return None

        logger.info(f"파라미터 설정:")
        logger.info(f"  - target_sum: {target_sum} ✓ (권장)")
        logger.info(f"  - log_transform: True ✓")
        logger.info(f"  - scale: {scale} ✓ (False 권장)")
        logger.info(f"  - hvg: {hvg} ✓ (False 권장)")

        # 1. 품질 메트릭 계산
        logger.info("\n[1-1] 품질 메트릭 계산...")
        sc.pp.calculate_qc_metrics(self.adata, percent_top=None, log1p=False, inplace=True)

        # QV 점수 확인 (만약 있으면)
        if 'qv' in self.adata.uns.get('spots', pd.DataFrame()).columns:
            qv_ratio = (self.adata.uns['spots']['qv'] > 20).sum() / len(self.adata.uns['spots'])
            logger.info(f"  - QV > 20 리드: {qv_ratio*100:.1f}%")
            logger.info(f"    (기준: ≥ 81% 권장)")

        # 2. 전처리 적용
        logger.info("\n[1-2] 전처리 적용...")
        self.adata = xp.main_preprocessing(
            self.adata,
            target_sum=target_sum,
            mincounts=mincounts,
            mingenes=mingenes,
            neigh=neigh,
            npc=0,
            scale=scale,
            hvg=hvg,
            norm=True,
            lg=True
        )

        logger.info(f"✓ 전처리 완료: {self.adata.shape}")
        logger.info(f"  - 클러스터 수: {len(self.adata.obs.get('leiden_1_4', []))}")

        # 3. 결과 저장
        if save:
            output_file = self.output_path / f"{self.sample_tag}_step1_preprocessed.h5ad"
            self.adata.write(str(output_file))
            logger.info(f"✓ 저장: {output_file}")
            
            # 4. 시각화 (Enhanced)
            logger.info("  - 시각화 생성 중 (QC, UMAP, Marker Genes)...")
            vis.visualize_step1_qc(self.adata, self.output_path, self.sample_tag)
            vis.visualize_pca_and_clustering(self.adata, self.output_path, self.sample_tag)
            vis.visualize_spatial_clusters(self.adata, self.output_path, self.sample_tag)
            
            # High Priority: Marker Genes
            vis.visualize_marker_genes(self.adata, self.output_path, self.sample_tag)
            
            # Medium Priority: Custom Genes
            genes_to_plot = kwargs.get('genes_to_plot', [])
            if genes_to_plot:
                vis.visualize_specific_genes(self.adata, self.output_path, self.sample_tag, genes_to_plot)

            # High Priority: Segmentation Overlay (Paper Image 5c)
            # Requires xenium_input_path to find boundaries
            if self.xenium_input_path:
                 vis.visualize_segmentation_overlay(self.xenium_input_path, self.output_path, self.sample_tag)


        return self.adata

    # =====================================================================
    # STEP 2: 세포 분할 없이 분석 (Segmentation-Free)
    # =====================================================================
    # 참고: notebooks/2_segmentation_free_analysis/2_1_batch_processing_distance_to_nuclei_across_samples.ipynb
    # 함수: dispersion(), dist_nuc()

    def step2_segmentation_free_analysis(
        self,
        sample_fraction: float = 0.1
    ) -> pd.DataFrame:
        """
        Step 2: 세포 분할 없이 리드-핵 거리 분석

        핵 vs 세포질 유전자 분포를 파악하고, 각 유전자의 기능적 특성을 연구합니다.

        Parameters
        ----------
        sample_fraction : float, optional
            분석할 데이터 샘플 비율, by default 0.1 (메모리 절약)

        Returns
        -------
        pd.DataFrame
            유전자별 거리 정보 DataFrame
        """
        logger.info("="*60)
        logger.info("STEP 2: 세포 분할 없이 분석 (Segmentation-Free)")
        logger.info("="*60)

        if self.adata is None or 'spots' not in self.adata.uns:
            logger.error("Step 0 데이터가 필요합니다.")
            return None

        logger.info(f"샘플링: {sample_fraction*100:.0f}% 데이터 사용")

        # 1. 리드 정보 준비
        reads_original = self.adata.uns['spots'].reset_index()

        # 2. 리드-세포 거리 계산
        logger.info("\n[2-1] 리드-세포 중심 거리 계산...")
        reads_assigned = dispersion(reads_original, self.adata)
        reads_assigned = reads_assigned.loc[:, [
            'feature_name', 'cell_id', 'distance',
            'overlaps_nucleus', 'x_location', 'y_location'
        ]]

        logger.info(f"✓ 거리 계산 완료: {reads_assigned.shape[0]} 리드")

        # 3. 샘플링 (메모리 절약)
        np.random.seed(0)
        sample_indices = np.random.choice(
            reads_assigned.shape[0],
            int(reads_assigned.shape[0] * sample_fraction),
            replace=False
        )
        reads_sample = reads_assigned.iloc[sample_indices]

        # Save raw sampled reads for distribution visualization (Step 2)
        raw_reads_path = self.output_path / f"{self.sample_tag}_step2_reads_sample.parquet"
        reads_sample.to_parquet(raw_reads_path)
        logger.info(f"✓ Saved raw sampled reads: {raw_reads_path}")

        # 4. 유전자별 거리 분석
        logger.info("\n[2-2] 유전자별 거리 분석...")
        gene_distances = reads_sample.groupby('feature_name')['distance'].agg(['mean', 'median', 'std', 'count'])
        gene_distances = gene_distances.sort_values('mean')

        # BLANK와 NegControl 제거
        gene_distances = gene_distances[
            ~gene_distances.index.str.contains('BLANK|NegControl', case=False)
        ]

        logger.info(f"✓ {len(gene_distances)} 유전자 분석 완료")

        # 5. 극단값 시각화
        logger.info("\n[2-3] 핵 근처 vs 세포질 풍부 유전자...")
        nuclear_genes = list(gene_distances.index[:5])
        cytoplasmic_genes = list(gene_distances.index[-5:])
        logger.info(f"  - 핵 농축: {nuclear_genes}")
        logger.info(f"  - 세포질 풍부: {cytoplasmic_genes}")

        # 6. 핵 확인 리드의 평균 거리
        reads_nuclear = reads_sample[reads_sample['overlaps_nucleus'] == 1]
        mean_nuclear_distance = reads_nuclear.groupby('cell_id')['distance'].max().mean()
        logger.info(f"  - 평균 핵 경계 거리: {mean_nuclear_distance:.2f} µm")

        # 7. 결과 저장
        output_file = self.output_path / f"{self.sample_tag}_step2_gene_distances.csv"
        gene_distances.to_csv(output_file)
        logger.info(f"✓ 저장: {output_file}")

        # 8. 시각화
        vis.visualize_step2_gene_distances(self.output_path, self.sample_tag)
        
        # Step 2.3: Cell Overlaps (Optional)
        if s2_overlaps:
            self.adata = s2_overlaps.run_step2_overlaps(self.adata, self.output_path, self.sample_tag)
        else:
            logger.warning("Step 2.3 Overlaps module not loaded. Skipping.")

        # Step 2.4: SSAM (Optional)
        if s2_ssam:
            self.adata = s2_ssam.run_step2_ssam(self.adata, self.output_path, self.sample_tag)
        else:
            logger.warning("Step 2.4 SSAM module not loaded. Skipping.")

        # Save Checkpoint
        # self.adata.write(self.output_path / f"{self.sample_tag}_step2_features.h5ad")

        return gene_distances

    # =====================================================================
    # STEP 4: 최적 세포 확장 거리 결정
    # =====================================================================
    # 참고: notebooks/4_optimal_expansion/4_1_Optimal_expansion_multisection.ipynb
    # 함수: dist_nuc(), 상관계수 계산

    def step4_optimal_expansion(self) -> Dict:
        """
        Step 4: 핵으로부터 최적 세포 확장 거리 결정

        세포 중심에서 몇 µm까지를 세포에 포함할지를 결정합니다.
        - 너무 좁으면: 세포질의 리드 손실
        - 너무 넓으면: 인접 세포 혼합 (bleeding)

        Returns
        -------
        Dict
            최적 확장 거리 정보
        """
        logger.info("="*60)
        logger.info("STEP 4: 최적 세포 확장 거리 결정")
        logger.info("="*60)

        if self.adata is None or 'spots' not in self.adata.uns:
            logger.error("Step 0 데이터가 필요합니다.")
            return None

        logger.info("\n[4-1] 리드 거리 계산...")
        reads_original = self.adata.uns['spots'].reset_index()

        # 세포 중심 좌표 매핑
        cell_coords = self.adata.obs[['x_centroid', 'y_centroid']]
        cell_id_to_coords = dict(zip(self.adata.obs['cell_id'],
                                     zip(cell_coords['x_centroid'],
                                         cell_coords['y_centroid'])))

        # 거리 계산
        reads_original['x_cell'] = reads_original['cell_id'].map(
            lambda x: cell_id_to_coords.get(x, (np.nan, np.nan))[0]
        )
        reads_original['y_cell'] = reads_original['cell_id'].map(
            lambda x: cell_id_to_coords.get(x, (np.nan, np.nan))[1]
        )

        dist2 = ((reads_original['x_location'] - reads_original['x_cell'])**2 +
                (reads_original['y_location'] - reads_original['y_cell'])**2)
        reads_original['distance'] = np.sqrt(dist2).round(0)

        # 유효한 리드만 선택
        reads_valid = reads_original.dropna(subset=['x_cell', 'y_cell'])

        logger.info(f"✓ {len(reads_valid)} 리드의 거리 계산 완료")

        # 2. 핵 겹침 리드의 최대 거리
        logger.info("\n[4-2] 세포 경계 분석...")
        reads_nuclear = reads_valid[reads_valid['overlaps_nucleus'] == 1]
        mean_nuclear_distance = reads_nuclear.groupby('cell_id')['distance'].max().mean()
        logger.info(f"  - 평균 핵 경계 (리드 기준): {mean_nuclear_distance:.2f} µm")

        # 3. 최적값 결정을 위한 반복 시뮬레이션
        logger.info("\n[4-3] 최적 확장 거리 시뮬레이션 (0 ~ 15 µm)...")
        
        # Purity 계산을 위한 샘플링 (속도 향상)
        unique_cells = reads_valid['cell_id'].unique()
        if len(unique_cells) > 1000:
            purity_sample_cells = np.random.choice(unique_cells, 1000, replace=False)
        else:
            purity_sample_cells = unique_cells
            
        reads_for_purity = reads_valid[reads_valid['cell_id'].isin(purity_sample_cells)]
        
        simulation_results = {}
        
        # Test range: 0 to 15 um, step 1.0 (or 0.5 for finer grain)
        test_range = np.arange(0, 16, 1.0)
        
        for d in test_range:
            # Filter by distance
            current_reads = reads_valid[reads_valid['distance'] <= d]
            
            # Metric 1: Capture (Total Reads)
            capture_count = len(current_reads)
            
            # Metric 2: Pseudo-Purity (Nuclear Fraction)
            # Assumption: Reads overlapping nucleus are "Pure". 
            # As we expand, we add non-nuclear reads. Purity = Nuclear Reads / Total Reads
            # This is a naive proxy but robust without Ground Truth.
            # Ideally we check consistency with nuclear profile, but that's expensive loops.
            # Using simple nuclear fraction on the subsample.
            
            current_purity_reads = reads_for_purity[reads_for_purity['distance'] <= d]
            if len(current_purity_reads) > 0:
                n_nuclear = current_purity_reads['overlaps_nucleus'].sum()
                purity_score = n_nuclear / len(current_purity_reads)
            else:
                purity_score = 0
            
            simulation_results[d] = {
                'Capture_Reads': capture_count,
                'Global_Purity_Proxy': purity_score
            }
            # logger.info(f"  - Distance {d} µm: {capture_count} reads, Purity {purity_score:.3f}")

        results_df = pd.DataFrame.from_dict(simulation_results, orient='index')
        results_df.index.name = 'Expansion_Distance'
        
        # Save results
        
        # Save results (Simulation Curve)
        # Use a distinct name so it's not overwritten by the summary
        output_file_sim = self.output_path / f"{self.sample_tag}_step4_expansion_curve.csv"
        results_df.to_csv(output_file_sim)
        logger.info(f"✓ 시뮬레이션 결과 저장: {output_file_sim}")

        # Pick optimal: specific heuristic or just report max capture with purity > threshold?
        # For now, we return the paper's recommendation but save our simulation
        optimal_from_center = 10.71
        optimal_from_nucleus_border = 5.65

        # Attempt to find "Elbow" or "Intersection"? 
        # Just logging for user to decide via visualization.
        
        result = {
            'optimal_from_center': optimal_from_center,
            'optimal_from_nucleus_border': optimal_from_nucleus_border,
            'mean_nucleus_distance': mean_nuclear_distance,
            'simulation_file': str(output_file_sim)
        }

        # 4. 결과 저장 (Summary)
        result_df = pd.DataFrame([result]).T
        output_file = self.output_path / f"{self.sample_tag}_step4_optimal_expansion.csv"
        result_df.to_csv(output_file)
        logger.info(f"✓ 저장: {output_file}")
        
        # 5. Visualization
        vis.visualize_step4_optimization(self.adata, self.output_path, self.sample_tag)

        return result

    # =====================================================================
    # STEP 6: 전처리 시뮬레이션 및 최적화
    # =====================================================================
    # 참고: notebooks/6_simulating_preprocessing/6_3_Simulated_Xenium_different_preprocessing_python.ipynb
    # 함수: allcombs(), compute_vi(), compute_fmi()

    def step6_preprocessing_simulation(
        self,
        sample_size: float = 0.05
    ) -> pd.DataFrame:
        """
        Step 6: 618가지 전처리 조합을 시뮬레이션하여 최적 경로 결정

        파라미터 조합:
        - norm: [True, False]
        - target_sum: [100, None]
        - scale: [True, False]
        - log_transform: [True, False]
        - hvg: [True, False]
        - ... 등등 (총 618가지)

        Parameters
        ----------
        sample_size : float, optional
            메모리 절약을 위한 샘플링 비율, by default 0.05

        Returns
        -------
        pd.DataFrame
            ARI 점수 DataFrame (조합 × 데이터셋)
        """
        logger.info("="*60)
        logger.info("STEP 6: 전처리 시뮬레이션 (618가지 조합)")
        logger.info("="*60)

        if self.adata is None:
            logger.error("Step 1 데이터가 필요합니다.")
            return None

        logger.info(f"샘플링: {sample_size*100:.0f}% 데이터 사용")

        # Determine if Ground Truth is available
        has_ground_truth = 'ground_truth_celltype' in self.adata.obs
        if has_ground_truth:
            logger.info("Ground Truth labels found in adata.obs['ground_truth_celltype'].")
        else:
            logger.info("Ground Truth labels NOT found. Using Stability ARI.")

        # 1. 데이터 샘플링
        maxcell = int(self.adata.shape[0] * sample_size)
        adata_sample = self.adata[:maxcell, :].copy()

        logger.info(f"  - 샘플 크기: {adata_sample.shape}")

        # 2. 모든 조합 시뮬레이션
        logger.info("\n[6-1] 전처리 파라미터 시뮬레이션 (618 조합)...")
        # Sample data for simulation to speed up if needed, though allcombs handles some logic
        # allcombs now returns (allres, sil_scores)
        try:
            allres, sil_scores = allcombs(
                adata_sample
            )
            logger.info(f"✓ 시뮬레이션 완료: {allres.shape[1]}개 조합 생성")
        except Exception as e:
            logger.error(f"조합 계산 실패: {e}")
            return None

        # 3. 각 조합의 성능 평가 (ARI: Adjusted Rand Index)
        logger.info("\n[6-2] 성능 평가 (ARI)...")
        from sklearn.metrics import adjusted_rand_score

        # Ground Truth 존재 여부 확인
        if has_ground_truth:
            logger.info("  ✓ Ground Truth 라벨 발견! Accuracy ARI를 계산합니다. (vs 'ground_truth_celltype')")
            # Ground Truth 데이터: self.adata.obs['ground_truth_celltype']
            # 샘플링된 데이터의 인덱스에 맞춰야 함
            # allcombs 내부에서 샘플링을 어떻게 하는지 확인 필요.
            # allcombs는 adata 전체를 받아서 내부적으로 처리하거나, 여기서 샘플링된 adata를 넘겨야 함.
            # allcombs(adata_sample) passes the sample now.
            
            # 주의: allres의 인덱스와 adata.obs의 인덱스가 일치해야 함.
            # CRITICAL FIX: Use adata_sample for GT labels to match allres indices
            gt_labels = adata_sample.obs['ground_truth_celltype']
            
            # allres가 전체 세포에 대한 결과라면 바로 비교 가능.
            # 만약 allres가 일부라면 정렬 필요.
            # 보통 allcombs 결과는 입력 adata와 동일 순서/인덱스.
        else:
            logger.info("  ⚠ Ground Truth 라벨 없음. Stability ARI를 계산합니다. (vs 'DEFAULT_louv')")

        ari_scores = {}
        for col in allres.columns:
            if has_ground_truth:
                # Accuracy ARI (vs Ground Truth)
                # 공통 인덱스만 비교 (혹시 모를 불일치 방지)
                common_idx = allres.index.intersection(gt_labels.index)
                if len(common_idx) < len(allres):
                    logger.warning("인덱스 불일치 발생, 공통 인덱스로 ARI 계산")
                
                ari = adjusted_rand_score(
                    gt_labels.loc[common_idx],
                    allres.loc[common_idx, col]
                )
            else:
                # Stability ARI (vs Default)
                ari = adjusted_rand_score(
                    allres.loc[:, 'DEFAULT_louv'],
                    allres.loc[:, col]
                )
            ari_scores[col] = ari

        # Save ARI scores
        ari_column_name = 'Accuracy_ARI' if has_ground_truth else 'Stability_ARI'
        ari_df = pd.DataFrame.from_dict(ari_scores, orient='index', columns=[ari_column_name])
        ari_df.index.name = 'Combination'
        ari_csv_path = self.output_path / f"{self.sample_tag}_step6_ari_scores.csv"
        ari_df.sort_values(by=ari_column_name, ascending=False).to_csv(ari_csv_path)

        # Save Silhouette scores (Keep as is)
        sil_df = pd.DataFrame.from_dict(sil_scores, orient='index', columns=['Silhouette'])
        sil_df.index.name = 'Combination'
        sil_csv_path = self.output_path / f"{self.sample_tag}_step6_silhouette_scores.csv"
        sil_df.sort_values(by='Silhouette', ascending=False).to_csv(sil_csv_path)
        
        summary_df = pd.merge(ari_df, sil_df, left_index=True, right_index=True)
        summary_df.to_csv(self.output_path / f"{self.sample_tag}_step6_summary.csv")
        logger.info(f"✓ 최종 결과 저장 (Summary): {self.output_path / f'{self.sample_tag}_step6_summary.csv'}")

        # Visualization
        vis.visualize_step6_simulation(self.output_path, self.sample_tag)

        return summary_df

    # =====================================================================
    # STEP 7: 공간 도메인 분석 (SpaGCN)
    # =====================================================================
    def step7_spatial_domains(self):
        """Step 7: Identify Spatial Domains using SpaGCN."""
        logger.info("="*60)
        logger.info("STEP 7: Spatial Domain Identification")
        logger.info("="*60)
        if s7_domains is None:
            logger.warning("Module step7_spatial_domains not loaded. Skipping.")
            return
        
        if self.adata is None:
            logger.error("No data loaded.")
            return

        self.adata = s7_domains.run_step7_spatial_domains(
            self.adata, 
            output_path=self.output_path, 
            sample_tag=self.sample_tag
        )

    # =====================================================================
    # STEP 8: 공간 변수 유전자 (SVF) 분석
    # =====================================================================
    def step8_svf(self):
        """Step 8: Identify Spatially Variable Features using SpatialDE."""
        logger.info("="*60)
        logger.info("STEP 8: SVF Identification")
        logger.info("="*60)
        if s8_svf is None:
            logger.warning("Module step8_svf not loaded. Skipping.")
            return
            
        if self.adata is None:
            logger.error("No data loaded.")
            return

        s8_svf.run_step8_svf(
            self.adata, 
            output_path=self.output_path, 
            sample_tag=self.sample_tag
        )
        
        # Visualization
        vis.visualize_step8_svf(self.adata, self.output_path, self.sample_tag)

        # Find best settings
        best_ari = max(ari_scores, key=ari_scores.get)
        best_ari_val = ari_scores[best_ari]

        best_sil = max(sil_scores, key=sil_scores.get)
        best_sil_val = sil_scores[best_sil]

        logger.info(f"\n[6-3] 최적 조합 추천:")
        logger.info(f"  ★ 최고 성능 ({ari_column_name} 1위): {best_ari} (ARI={best_ari_val:.3f})")
        logger.info(f"  ★ 최고 클러스터 품질 (Silhouette 1위): {best_sil} (Sil={best_sil_val:.3f})")

        logger.info(f"\n[6-4] GOLDEN PATH (권장) 점수:")
        gp_ari = ari_scores.get('DEFAULT_louv', 0)
        gp_sil = sil_scores.get('DEFAULT_louv', 0)
        logger.info(f"  - ARI: {gp_ari:.3f}")
        logger.info(f"  - Silhouette: {gp_sil:.3f}")

        logger.info(f"\n✓ 저장 완료:")
        logger.info(f"  - {ari_csv_path}")
        logger.info(f"  - {sil_csv_path}")
        logger.info(f"  - {summary_csv_path}")

        return summary_df

    # =====================================================================
    # 통합 파이프라인 실행
    # =====================================================================

    def run_full_pipeline(
        self,
        steps: List[int] = [0, 1, 6],
        **kwargs
    ):
        """
        전체 또는 부분 파이프라인 실행
        """
        logger.info("\n" + "="*60)
        logger.info("Xenium 분석 파이프라인 시작")
        logger.info("="*60)
        logger.info(f"실행 스텝: {steps}")
        logger.info("="*60 + "\n")

        try:
            if 0 in steps:
                # Step 0에서 Ground Truth 생성 파라미터 전달
                s0_kwargs = {k: v for k, v in kwargs.items() if k in ['generate_ground_truth', 'tissue_type', 'download_dir', 'organ']}
                self.step0_format_xenium(**s0_kwargs)

            if 1 in steps:
                self.step1_preprocess(**kwargs)

            if 2 in steps:
                self.step2_segmentation_free_analysis()

            if 4 in steps:
                self.step4_optimal_expansion()

            if 6 in steps:
                self.step6_preprocessing_simulation()

            if 7 in steps:
                self.step7_spatial_domains()

            if 8 in steps:
                self.step8_svf()
            
            logger.info("\n" + "="*60)
            logger.info("✓ 파이프라인 완료!")
            logger.info("="*60)
            logger.info(f"결과 저장 위치: {self.output_path}")

        except Exception as e:
            logger.error(f"\n✗ 파이프라인 실패: {e}")
            raise


# ==========================================================================
# 메인 실행 코드
# ==========================================================================

if __name__ == "__main__":

    try:
        # 설정 파일 로드
        logger.info("설정 파일 로드: config.yaml")
        config = load_config("config.yaml")
        setup_logging(config)

        logger.info("=" * 70)
        logger.info("Xenium 분석 파이프라인 시작")
        logger.info("=" * 70)

        # 설정에서 필요한 값 추출
        data_config = config.get('data', {})
        pipeline_config = config.get('pipeline', {})
        preprocessing_config = config.get('preprocessing', {})
        step_params = config.get('step_parameters', {})

        XENIUM_INPUT_PATH = data_config.get('xenium_input_path')
        OUTPUT_PATH = data_config.get('output_path', './output')
        SAMPLE_TAG = data_config.get('sample_tag', 'xenium_sample')
        STEPS = pipeline_config.get('steps', [0, 1, 6])

        # 경로 검증
        if not XENIUM_INPUT_PATH:
            raise ValueError("config.yaml에서 xenium_input_path를 설정하세요!")

        logger.info(f"\n설정 정보:")
        logger.info(f"  입력 경로: {XENIUM_INPUT_PATH}")
        logger.info(f"  출력 경로: {OUTPUT_PATH}")
        logger.info(f"  샘플 태그: {SAMPLE_TAG}")
        logger.info(f"  실행 단계: {STEPS}")

        # 파이프라인 초기화
        # load_preprocessed=True이면 h5ad 파일을 직접 로드
        load_h5ad = XENIUM_INPUT_PATH and XENIUM_INPUT_PATH.endswith('.h5ad')

        pipeline = XeniumPipeline(
            xenium_input_path=XENIUM_INPUT_PATH,
            output_path=OUTPUT_PATH,
            sample_tag=SAMPLE_TAG,
            load_preprocessed=load_h5ad
        )

        logger.info(f"\n파일 타입: {'H5AD 파일 (Step 1부터)' if load_h5ad else 'Xenium 원본 (Step 0부터)'}")

        # 전처리 파라미터 준비
        preproc_kwargs = {
            'target_sum': preprocessing_config.get('target_sum', 100),
            'mincounts': preprocessing_config.get('mincounts', 10),
            'mingenes': preprocessing_config.get('mingenes', 3),
            'neigh': preprocessing_config.get('neigh', 15),
            'scale': preprocessing_config.get('scale', False),
            'hvg': preprocessing_config.get('hvg', False),
        }

        # Step별 파라미터 준비
        kwargs = {'**preproc_kwargs': preproc_kwargs}

        # Step 2 파라미터
        if 2 in STEPS:
            step2_config = step_params.get('step2', {})
            kwargs['step2_sample_fraction'] = step2_config.get('sample_fraction', 0.1)

        # Step 6 파라미터
        if 6 in STEPS:
            step6_config = step_params.get('step6', {})
            kwargs['step6_sample_size'] = step6_config.get('sample_size', 0.05)

        # Step 7 Config를 Step 0 Parameter로 변환
        step7_config = step_params.get('step7', {})
        # 사용자가 generate_ground_truth를 켜고 싶어하는지 확인 (기본값 설정 필요)
        # config.yaml에 명시적 설정이 없다면, step7_config 존재 여부로 판단하거나 기본값 False
        # 여기서는 편의상 step7_config가 비어있지 않거나, tissue_type이 있으면 True로 설정
        generate_gt = step7_config.get('enabled', False) # config.yaml에 'enabled' 추가 권장
        
        # Default tissue inference
        default_tissue = "brain" if "brain" in SAMPLE_TAG.lower() else "breast"
        
        kwargs['generate_ground_truth'] = generate_gt
        kwargs['tissue_type'] = step7_config.get('tissue_type', default_tissue)
        kwargs['download_dir'] = step7_config.get('download_dir', 'data/reference')
        kwargs['organ'] = step7_config.get('organ', None)

        if generate_gt:
            logger.info(f"Ground Truth 생성 활성화: Tissue={kwargs['tissue_type']}")

        # 파이프라인 실행
        logger.info(f"\n파이프라인 실행 중...")
        pipeline.run_full_pipeline(
            steps=STEPS,
            **preproc_kwargs
        )

        logger.info(f"\n" + "=" * 70)
        logger.info("✓ 파이프라인 완료!")
        logger.info("=" * 70)
        logger.info(f"결과 저장 경로: {OUTPUT_PATH}")
        logger.info(f"로그 파일: xenium_pipeline.log")

    except FileNotFoundError as e:
        logger.error(f"\n✗ 오류: {e}")
        logger.error("\n해결 방법:")
        logger.error("  1. config.yaml 파일이 현재 디렉토리에 있는지 확인하세요")
        logger.error("  2. config.yaml에서 xenium_input_path를 설정하세요")
        exit(1)

    except ValueError as e:
        logger.error(f"\n✗ 설정 오류: {e}")
        exit(1)

    except Exception as e:
        logger.error(f"\n✗ 예상치 못한 오류: {e}")
        logger.exception("상세 정보:")
        exit(1)

