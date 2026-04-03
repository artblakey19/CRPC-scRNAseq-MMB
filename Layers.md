# Seurat v5 IntegrateLayers 워크플로우 정리

## 1. Layer 구조란?

Seurat v5에서는 RNA assay 내부에 **환자(sample)별로 분리된 layer**가 존재한다.
예: `counts.CRPC1`, `counts.CRPC2`, `counts.CRPC3`

- `split()` → layer 분리
- `JoinLayers()` → layer 합침

---

## 2. IntegrateLayers가 하는 일

- 분리된 layer들을 **서로 비교**하여 CCA anchor를 찾고, batch effect를 보정한 저차원 embedding(`integrated.dr`)을 생성한다.
- **원본 count matrix 자체를 수정하지 않는다.** 보정된 embedding만 만들어줌.

### 따라서 반드시 복수 layer가 필요

| 함수 | 입력 | layer 분리 필요 여부 |
|---|---|---|
| SCTransform | 셀 × gene count matrix | 불필요 (셀 단위 연산) |
| RunPCA | SCT residual matrix | 불필요 (matrix 단위 연산) |
| **IntegrateLayers** | **layer 간 비교** | **필수 (복수 layer)** |

---

## 3. Subset 후 Re-integration 시 주의사항

전체 세포에서 `JoinLayers()`를 이미 수행한 뒤 `subset()`하면 **단일 layer** 상태가 된다.
이 상태에서 `IntegrateLayers()`를 호출하면 비교 대상이 없어 에러 발생:

```
Error in UseMethod(generic = "Assays", object = object) :
  no applicable method for 'Assays' applied to an object of class "NULL"
```

### 해결: subset 후 split으로 layer 재분리

```r
epi <- subset(combined_CRPC, subset = celltype == "Epithelial")
epi[["RNA"]] <- split(epi[["RNA"]], f = epi$orig.ident)  # ← 이 줄이 핵심

epi <- SCTransform(epi, verbose = FALSE)
epi <- RunPCA(epi, verbose = FALSE)
epi <- IntegrateLayers(object = epi, method = CCAIntegration,
                        normalization.method = "SCT", verbose = FALSE)
epi[["RNA"]] <- JoinLayers(epi[["RNA"]])
```

---

## 4. Re-integration을 하는 이유

전체 세포 대상으로 수행한 CCA integration은 **모든 세포 타입의 변이를 포함한 상태에서** 환자 간 공통 변동축을 찾은 것이다.

Epithelial만 subset한 뒤 re-integration하면:
- CCA가 **epithelial 내부의 변이만을 대상으로** 환자 간 공통축을 새로 찾음
- Epithelial 내부의 미세한 subtype 차이 (luminal, basal, neuroendocrine 등)가 더 잘 드러남
- Batch correction이 해당 세포군에 맞게 더 정밀해짐

---

## 5. JoinLayers가 Integration 뒤에 필요한 이유

IntegrateLayers 후에도 RNA assay의 count/data layer는 **여전히 환자별로 분리된 상태**이다.

| Downstream 분석 | 사용하는 것 | JoinLayers 없이 가능? |
|---|---|---|
| FindNeighbors / FindClusters | `integrated.dr` (embedding) | ✅ 가능 |
| RunUMAP | `integrated.dr` (embedding) | ✅ 가능 |
| **FindMarkers / FindAllMarkers** | **원본 count/data matrix** | ❌ 불가 → 에러 |

→ `JoinLayers()`는 분리된 expression matrix를 하나로 합쳐서 DEG 분석을 가능하게 해주는 단계.