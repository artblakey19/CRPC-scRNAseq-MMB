# copyKAT 아카이브 (2026-07-24)

CNV 추정을 copyKAT에서 **철수**하고 아카이브. 추후 **Numbat**으로 재시도 예정.

## 내용물

| 경로 | 설명 |
|---|---|
| `scripts/08_CopyKAT/01_Run_CopyKAT.R` | per-sample copyKAT 실행 (마지막 상태 = 면역-only anchor) |
| `scripts/08_CopyKAT/02_CopyKAT_UMAP_Epithelial.R` | 예측을 상피 UMAP(stage 05)에 투사 |
| `Results/08_CopyKAT/` | 마지막 실행(면역-only anchor) 산출물 18GB |
| `logs/08_copykat_CRPC{1,2,3}.log` | 실행 로그 |

⚠️ `Results/08_CopyKAT/`는 **면역-only anchor 실행 결과**로, 그 이전의 raw-counts canonical(비상피 전체 anchor) 산출물을 덮어쓴 상태다. canonical 산출물 백업은 없다(`Results/`는 gitignore). copyKAT은 `CNA.MCMC`가 `MCMCpack::MCpoissongamma(mc=1000)` 몬테카를로 표본에 KS test를 걸고 그게 `mclapply` 안에서 도는데 `RNGkind("L'Ecuyer-CMRG")`가 없어 자식 seed가 시간+PID로 잡히므로, 래퍼의 `set.seed(42)`로도 재실행 시 동일 수치 재현 불가.

**소실된 canonical(비상피 전체 anchor, raw counts) 수치 — 기록 보존용**:

| | CRPC1 | CRPC2 | CRPC3 |
|---|---|---|---|
| aneuploid | **0** | 1,314 | 5,189 |

aneuploid는 상피에만 나왔고 비상피는 전부 diploid였다(anchor 설계대로의 sanity check 통과). ※ 이 산출물은 `n.cores=16` 시절 생성분이고 현행 스크립트는 30이라, 되돌려 재실행해도 근사조차 보장되지 않는다.

## 철수 사유 — CRPC1 baseline 실패

마지막 실행에서 anchor를 비상피 전체 → 면역세포(T/NK + Phagocytes + Mast)로 좁힌 결과, CRPC1에서 판정이 붕괴.

aneuploid 비율(celltype × sample):

| celltype | CRPC1 | CRPC2 | CRPC3 |
|---|---|---|---|
| Smooth muscle | **37.2%** | 0% | 0% |
| Endothelial | **17.5%** | 0.1% | 0% |
| Fibroblast | **5.7%** | 0% | 0% |
| anchor 면역세포 자체 | **9.1%** | 0% | 0% |
| Epithelial | 68.3% | 22.0% | 27.6% |
| `not.defined` (전체) | 881 (7.3%) | 235 (2.5%) | 84 (1.0%) |

**내부 모순**: diploid baseline으로 직접 투입한 anchor 면역세포 157개(9.1%)가 aneuploid로 판정됨. 정상 기질(평활근·내피)도 대량 aneuploid. CRPC2/3는 정확히 0%로 깨끗 → CRPC1에만 기술적 위양성 바닥값이 존재하며, 상피 68.3%도 같은 편향에 부풀려진 값이라 환자 간 비교 불가.

**동일 클러스터의 환자별 불일치** (aneuploid %):

| annotation | CRPC1 | CRPC2 | CRPC3 |
|---|---|---|---|
| BE 4 | 78.8% | 0.5% | 6.8% |
| BE 6 | 72.1% | 0.9% | 0% |
| OE | 74.0% | 0.3% | 2.4% |

동일 전사체 클러스터가 환자별로 정반대 → per-sample baseline drift.

**원인** (anchor 수 부족은 아님 — 로그상 baseline 사용 normal은 CRPC1 917개로 최다, CRPC2 362, CRPC3 244):

copyKAT은 `cl.aneuploid <- which(cl.ID == min(cl.ID))`로 **"종양이 있는가"를 묻지 않고 무조건 종양 극을 지명**한다. 비상피 전체를 anchor로 쓰면 모든 정상 클러스터가 anchor를 포함해 `cl.ID`가 "알려진 정상 vs 정체불명 상피"를 가르지만, 면역세포만 쓰면 정상 기질·정상 상피가 전부 `cl.ID≈0`이 되어 축이 노이즈화되고 극 지명이 임의가 된다. CRPC1은 종양 함량이 낮을 가능성이 있어(biopsy histology 거의 benign) 이 실패에 특히 취약하다.

이는 2026-07-23 **T/NK-only anchor 실험의 재현**이다. 당시 수치: CRPC1 상피 aneuploid 0→5,528(70%), anchor T/NK 자신 9.8% aneuploid, smooth muscle 37%·endothelial 17%. 이번(면역 3종 anchor) 68.3% / 9.1% / 37.2% / 17.5% — 사실상 동일. 그때 내린 **"좁은 anchor 재시도 금지"** 결론이 그대로 재확인됐다.

## Numbat 재시도 시 유의

- **좁은 anchor 금지.** 면역-only / T/NK-only / ionocyte-anchor 모두 실패 이력.
- **CRPC1이 진짜 시험대.** 종양 함량이 낮을 가능성이 있어, 기질·면역세포가 정상(diploid)으로 나오는지를 먼저 sanity check로 볼 것. 상피 aneuploid 비율보다 이 지표가 신뢰도를 가른다.
- Numbat은 haplotype-aware(allele + expression)라 baseline 강제 이분할 문제에 덜 취약할 것으로 기대되나, phasing 참조 패널 준비가 필요.
- 입력은 raw counts canonical 유지(decontX 입력판은 aneuploid를 부풀림: CRPC2 3945 → raw 1314).

## 변경되지 않는 규칙

CNV 결과는 malignancy / subclone / trajectory 어떤 argument에도 근거로 사용하지 않는다(supporting evidence로도 금지). 위 수치는 **진단용**으로만 기록.
