# PCa 종양-특이 전사체 마커 아카이브 (2026-07-24)

copyKAT 철수([[../copykat_2026-07-24]]) 후 CNV 대체로 시도한 **전립선암 특이 전사체 FeaturePlot**. 사용자 결정으로 아카이브하고 **Numbat**으로 넘어감.

## 내용물

| 경로 | 설명 |
|---|---|
| `scripts/08_Malignancy/01_PCa_Tumor_Transcripts_FeaturePlot.R` | 유전자당 4패널(annotation + CRPC1/2/3 split) FeaturePlot, assay=SCT, viridis, 3환자 공통 색스케일/프레임 |
| `Results/08_Malignancy/Figures/` | 9개 유전자 PNG |

## 아카이브 사유 — DNPC 코호트에 부적합

플롯한 마커(PCA3·AMACR·SCHLAP1·ERG·TMEFF2·GOLM1·MYC·NKX3-1·PTEN)가 **전부 AR축/luminal 선암 마커**라, AR-negative가 정의인 DNPC 코호트에서는 LE(ARPC) 잔존물에만 켜지고 본체 상피에선 대부분 꺼짐. malignancy 판별에 쓸 신호가 약함.

발현 세포 비율(전체 19,930 상피):

| 유전자 | 양성세포 | % | 비고 |
|---|---|---|---|
| PCA3 | 45 | 0.23% | 사실상 백지, 전부 CRPC2 |
| TMEFF2 | 60 | 0.30% | 사실상 백지, 전부 CRPC1 |
| ERG | 748 | 3.75% | 희박, CRPC3 편중 |
| AMACR | 972 | 4.88% | 확산적, 두드러진 클러스터 없음 |
| SCHLAP1 | 1,464 | 7.35% | Hillock1 30.7% |
| GOLM1 | 2,799 | 14.0% | LE(ARPC) 48.6% (유일하게 뚜렷) |
| NKX3-1 | 4,385 | 22.0% | 소실 마커 — dropout과 구분 불가 |
| PTEN | 7,698 | 38.6% | 소실 마커 — 변별력 낮음 |
| MYC | 8,904 | 44.7% | 광범위, 변별력 낮음 |

추가 한계: NKX3-1/PTEN 같은 소실 마커는 10x 3' dropout과 생물학적 소실을 세포 단위로 구분 불가. AR-driven signature가 DNPC에서 ARPC 외 음수인 기존 결론([[project_pca_modulescore_dnpc_limitation]])과 같은 맥락.

## 남는 방향 (미실행)

DNPC 본체에는 양성 종양 마커가 없음(정의상 AR-null + NE-null 소거법). 논의만 하고 안 만든 후속:
1. NE 소거 확증 패널(CHGA/CHGB/SYP/INSM1/NCAM1/ASCL1) — double-negative의 N쪽 검증
2. 재현되는 양성 프로그램(basal-squamous 케라틴 + CXCL8 염증) — DNPC state
필요 시 여기 스크립트를 참고해 재작성 가능.

## 다음 단계

CNV는 **Numbat**(haplotype-aware, allele+expression)으로 재시도. copyKAT의 baseline 강제 이분할 문제에 덜 취약할 것으로 기대. CRPC1의 기질·면역이 diploid로 나오는지를 1차 sanity check로 볼 것.
