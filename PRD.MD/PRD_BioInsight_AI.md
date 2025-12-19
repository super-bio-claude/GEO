# Product Requirements Document (PRD)
# BioInsight AI

**Version:** 1.0  
**Date:** 2024.12.19  
**Author:** Product Team  
**Status:** Draft

---

## 1. Executive Summary

BioInsight AI는 바이오·헬스케어 연구자들을 위한 AI 기반 웹 연구 지원 플랫폼입니다. 복잡한 데이터 분석과 분산된 문헌 해석 과정을 통합하여, 논문 분석, RNA-seq 데이터 분석, 머신러닝 기반 예측 모델을 하나의 일관된 분석 흐름으로 제공합니다.

### 제품 개요

| 항목 | 내용 |
|------|------|
| **타겟 사용자** | 대학/연구소 연구원, 바이오 스타트업, 제약사 R&D |
| **핵심 가치** | 데이터 분석 + 논문 지식 + ML 예측의 원스톱 통합 |
| **차별화 포인트** | 초보자도 전문가 수준의 분석 가능, LLM 기반 해석 |

---

## 2. Problem Statement

### 2.1 연구자 Pain Points

| 문제 영역 | 구체적 상황 | 영향 |
|----------|------------|------|
| **논문 해석 부담** | 복잡한 구조, 방대한 내용, 전문 용어로 1편당 2-4시간 소요 | 연구 시간 50% 이상 문헌 리뷰에 소모 |
| **RNA-seq 분석 장벽** | DESeq2, apeglm, clusterProfiler 등 전문 도구 러닝커브 | Wet Lab 연구자의 데이터 활용 제한 |
| **정보 통합 부재** | 논문, DEG 결과, pathway를 연결하여 해석하는 도구 없음 | 단편적 분석, 인사이트 도출 어려움 |
| **정보 과부하** | PubMed 기준 매일 수백 편의 신규 논문 발표 | 최신 연구 트렌드 파악 불가 |

### 2.2 기존 솔루션의 한계

| 기존 도구 | 제공 기능 | 미비점 |
|----------|----------|--------|
| DESeq2 | Differential Expression 통계 분석 | ML 없음, 코딩 필수 |
| GEO2R | 웹 기반 DE 분석 | ML 없음, 논문 해석 없음 |
| Galaxy | 워크플로우 기반 NGS 분석 | 논문 지식 통합 없음 |
| Tempus/Owkin | AI 기반 정밀의료 분석 | 일반 연구자 접근 불가 |

---

## 3. Solution Overview

### 3.1 제품 정의

BioInsight AI는 연구자가 논문(PDF) 또는 RNA-seq 데이터를 업로드하면 AI가 자동으로 분석하고, 핵심 정보를 추출하며, 관련 연구까지 추천하는 웹 기반 분석 플랫폼입니다. PubMed 논문 RAG 기능을 통합하여 "실험 결과 + 최신 논문 근거"를 동시에 제시합니다.

### 3.2 핵심 기능 구성

#### 3.2.1 논문 분석 모듈

- PDF 업로드 → 자동 파싱 (제목/초록/본문 구조 자동 구분)
- BioBERT/PubMedBERT 기반 문헌 임베딩 생성
- LLM 기반 핵심 요약 및 Q&A 제공
- 유사 논문 추천 (ChromaDB 벡터 검색)

#### 3.2.2 RNA-seq 분석 자동화

- Count Matrix & Metadata 업로드
- DESeq2 기반 통계 분석 자동 수행 (R→Python 자동 번역)
- 시각화: Volcano Plot, PCA, Heatmap 자동 생성
- apeglm/ashr 기반 로그폴드 조정
- Pathway/Gene/Ontology 분석 (clusterProfiler 연동)

#### 3.2.3 ML 예측 모델

- 선택된 유전자 Signature로 ML 모델 자동 학습
- 진단/그룹 분류 예측 (XGBoost, Random Forest)
- ROC-AUC, Precision/Recall 검증
- SHAP 기반 Feature Importance 시각화

#### 3.2.4 AI 연구 어시스턴트

- 논문 내용 + DEG 결과 + Pathway 정보를 종합하여 자동 답변하는 통합 챗봇
- "이 논문 결과와 내 RNA-seq 데이터는 어떤 관계인가?"와 같은 복합 질문 처리
- 추가 실험 아이디어 제안 기능

---

## 4. User Flow

### 4.1 논문 분석 Flow

```
PDF 업로드 → 자동 파싱 및 구조화 → BioBERT 임베딩 생성 → 핵심 요약 및 Q&A 제공 → 유사 논문 추천
```

### 4.2 RNA-seq 분석 Flow

```
1. Count Matrix 및 실험 Metadata 업로드
2. DESeq2 기반 자동 분석 수행 (정규화, Batch Correction)
3. 시각화 제공 (Volcano, PCA, Heatmap)
4. DEG 기반 Pathway 및 Gene Ontology 분석 결과 생성
5. ML 기반 예측 모델 실행 (샘플 특성/조건에 대한 예측 수행)
6. 해석 중심의 분석 리포트 제공
7. 예측값 + 해석 결과 리포트 제공
```

### 4.3 통합 분석 Flow

AI 연구 어시스턴트는 RNA-seq 기반 DEG 및 pathway 분석 결과를 입력으로 받아, 관련 문헌의 유전자·기전·표현형 정보를 정렬·비교하여 **"내 실험 데이터가 기존 연구를 재현하는지, 확장하는지, 혹은 새로운 관점을 제시하는지"를 자동으로 분석·요약**합니다.

---

## 5. Technical Architecture

### 5.1 기술 스택

| 레이어 | 기술 |
|--------|------|
| **Backend** | Python, FastAPI, rpy2 (DESeq2 호출), sklearn, XGBoost |
| **Frontend** | Streamlit (MVP) → React + Tailwind (정식 버전) |
| **Database** | PostgreSQL (결과 저장), ChromaDB (문헌 임베딩) |
| **AI/ML** | PubMedBERT/BioBERT, GPT-4o/Claude/Gemini, LangChain (RAG) |
| **인프라** | AWS/GCP, Docker, CI/CD Pipeline |

### 5.2 RNA-seq 분석 파이프라인

실제 질병 예측 모델 구축을 위한 표준 파이프라인:

1. **데이터 확보:** RNA-seq Count Matrix (n × genes 형태)
2. **정규화:** TPM, FPKM, VST (Variance Stabilizing Transformation)
3. **Batch Correction:** ComBat, Limma removeBatchEffect
4. **DE 분석:** DESeq2/edgeR로 질병군 vs 정상군 유전자 차이 확인
5. **Feature Engineering:** 유전자 Signature 생성, PCA/Autoencoder 활용
6. **ML 모델 학습:** 진단/중증도/클래스 분류
7. **검증:** ROC-AUC, Precision/Recall, SHAP Feature Importance

### 5.3 시스템 아키텍처 다이어그램

```
┌─────────────────────────────────────────────────────────────────────┐
│                           Frontend Layer                             │
│                    Streamlit (MVP) / React + Tailwind                │
└─────────────────────────────────┬───────────────────────────────────┘
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────────┐
│                           API Gateway                                │
│                             FastAPI                                  │
└────────┬─────────────────────┬──────────────────────┬───────────────┘
         │                     │                      │
         ▼                     ▼                      ▼
┌─────────────────┐  ┌─────────────────┐  ┌─────────────────────────┐
│  논문 분석 모듈  │  │ RNA-seq 분석   │  │    ML 예측 모듈        │
│  - PDF Parser   │  │ - DESeq2(rpy2) │  │    - XGBoost           │
│  - BioBERT      │  │ - Visualization│  │    - Random Forest     │
│  - LangChain    │  │ - clusterProfiler│ │    - SHAP              │
└────────┬────────┘  └────────┬────────┘  └───────────┬─────────────┘
         │                     │                      │
         ▼                     ▼                      ▼
┌─────────────────────────────────────────────────────────────────────┐
│                         Data Layer                                   │
│         PostgreSQL (Results)  │  ChromaDB (Embeddings)              │
└─────────────────────────────────────────────────────────────────────┘
```

---

## 6. Differentiation Strategy

### 6.1 핵심 차별점

| 차별화 요소 | 세부 내용 |
|------------|----------|
| **LLM 기반 논문 해석 자동화** | 복잡한 학술 논문을 자동 요약하고 Q&A 제공 |
| **ML 모델 자동 생성** | 코딩 없이 유전자 데이터로 예측 모델 생성 |
| **실험 + 문헌 통합 리포트** | 분석 결과와 관련 논문 근거를 하나의 보고서로 제공 |
| **민주화된 접근성** | Wet Lab 연구자도 전문가 수준의 분석 가능 |

### 6.2 핵심 가치 제안

> *"단순 분석 툴이 아니라, 데이터 분석 + 논문 지식 통합 + ML 실험까지 일괄 제공하는 최초의 바이오 AI 플랫폼"*

---

## 7. Expected Outcomes

### 7.1 정량적 효과

| 지표 | 기존 방식 | BioInsight AI |
|------|----------|---------------|
| 논문 1편 분석 시간 | 2-4시간 | **10-30분** |
| RNA-seq 분석 소요 기간 | 1-2주 (외부 의뢰 시) | **수 시간 내** |
| ML 모델 구축 | 별도 전문가 필요 | **자동 생성** |

### 7.2 정성적 효과

- **연구 생산성 향상:** 데이터 해석 시간 단축으로 핵심 연구에 집중
- **디지털 전환 가속:** Dry Lab 역량 강화 및 실험-분석 연계 효율화
- **연구 품질 향상:** 문헌 근거 기반의 체계적인 분석 및 해석
- **AI 연구 조력자:** 연구자 개인 맞춤형 AI 어시스턴트로 발전 가능

---

## 8. Development Roadmap

| 단계 | 기간 | 목표 | 주요 산출물 |
|------|------|------|------------|
| **Phase 1** | 1-3개월 | MVP 개발 | Streamlit 기반 핵심 기능 프로토타입 |
| **Phase 2** | 4-6개월 | Beta 출시 | React UI, 전체 파이프라인 통합 |
| **Phase 3** | 7-12개월 | 정식 출시 | 고급 기능, 엔터프라이즈 버전 |

### Phase 1 세부 마일스톤 (MVP)

- [ ] 논문 PDF 업로드 및 파싱 기능
- [ ] 기본 요약 및 Q&A 기능
- [ ] RNA-seq Count Matrix 업로드
- [ ] DESeq2 기반 분석 파이프라인
- [ ] 기본 시각화 (Volcano, PCA)

### Phase 2 세부 마일스톤 (Beta)

- [ ] React 기반 UI 전환
- [ ] ChromaDB 기반 유사 논문 추천
- [ ] Pathway/GO 분석 통합
- [ ] ML 예측 모델 모듈
- [ ] 통합 리포트 생성

### Phase 3 세부 마일스톤 (정식 출시)

- [ ] 엔터프라이즈 보안 기능
- [ ] 팀 협업 기능
- [ ] API 제공
- [ ] 커스텀 모델 학습 지원

---

## 9. Success Metrics

### 9.1 핵심 KPI

| 지표 | 목표 (Phase 1) | 목표 (Phase 3) |
|------|---------------|---------------|
| MAU (Monthly Active Users) | 100명 | 5,000명 |
| 분석 완료율 | 70% | 90% |
| 평균 분석 시간 | 30분 이내 | 15분 이내 |
| 사용자 만족도 (NPS) | 30+ | 50+ |
| 재사용률 | 40% | 70% |

---

## 10. Risks & Mitigations

| 리스크 | 영향도 | 발생 가능성 | 대응 방안 |
|--------|--------|------------|----------|
| RNA-seq 분석 정확도 이슈 | 높음 | 중간 | 전문가 검증, 벤치마크 데이터셋 활용 |
| LLM 환각(Hallucination) | 높음 | 중간 | RAG 기반 근거 제시, 출처 명시 |
| 사용자 학습 곡선 | 중간 | 높음 | 온보딩 튜토리얼, 샘플 데이터 제공 |
| 인프라 비용 증가 | 중간 | 중간 | 캐싱 전략, 사용량 기반 과금 |

---

## 11. Conclusion

BioInsight AI는 기존에 분산되어 제공되던 RNA-seq 분석, 기능 해석, 문헌 탐색, 예측 분석 도구들을 하나의 플랫폼으로 통합하고, LLM 기반의 지능형 해석 기능과 머신러닝 기반 예측 모델을 결합함으로써 **연구자가 전문가의 분석 사고 과정을 보다 쉽게 재현할 수 있는 바이오 데이터 분석 환경을 제공합니다.**

데이터 중심의 연구가 필수가 된 현대 바이오 연구 환경에서, BioInsight AI는 반복적인 분석과 문헌 탐색에 소요되는 시간을 줄이고, 실험 결과 해석의 일관성과 깊이를 향상시켜 연구 품질을 높이는 핵심 연구 지원 도구로 활용될 것입니다.

---

## Appendix

### A. 용어 정의

| 용어 | 정의 |
|------|------|
| RNA-seq | RNA 시퀀싱, 전사체 분석 기술 |
| DEG | Differentially Expressed Genes, 차등 발현 유전자 |
| DESeq2 | RNA-seq 데이터 분석을 위한 R 패키지 |
| RAG | Retrieval-Augmented Generation, 검색 증강 생성 |
| Pathway | 생물학적 경로, 유전자들의 상호작용 네트워크 |
| Gene Ontology | 유전자 기능 분류 체계 |
| SHAP | SHapley Additive exPlanations, ML 모델 해석 기법 |

### B. 참고 문헌

- DESeq2: Love, M.I., Huber, W., Anders, S. (2014)
- BioBERT: Lee, J., et al. (2020)
- LangChain Documentation

---

*Document Version History*

| 버전 | 날짜 | 변경 내용 | 작성자 |
|------|------|----------|--------|
| 1.0 | 2024.12.19 | 최초 작성 | Product Team |
