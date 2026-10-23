# Validación clínica de un fármaco reposicionado en una enfermedad ultra-rara (<50 pacientes): evidencia, diseños, vías regulatorias y precedentes

Nota de método: ~17 llamadas de búsqueda/fetch. Los snippets de búsqueda se citan con la URL de origen; lo marcado "[de memoria, sin verificar]" no se pudo confirmar en esta sesión y debe comprobarse antes de citarlo en el informe. Nada de esto contiene datos del paciente.

## 1. Escalera de evidencia (células de paciente, organoides/iPSC, modelos de ratón existentes, primer uso en paciente)

### Takeaway
No existe una escalera normada única: en ultra-rara se usa en la práctica una cadena "mecanismo genético conocido -> prueba en células del propio paciente -> (opcional) modelo animal -> uso en el paciente con biomarcador farmacodinámico". El cribado fenotípico de fármacos aprobados sobre células primarias del paciente es la estrategia de reposicionamiento más directa, y la FDA acepta explícitamente datos de modelos no clínicos y farmacodinámicos como evidencia confirmatoria.

### Cited Findings
- Revisión: el cribado fenotípico + reposicionamiento con células primarias de paciente y fenotipos celulares que explican los síntomas es descrito como la estrategia más directa y relevante en trastornos ultra-raros — [Drug repurposing and phenotypic screening: innovative strategies for treating ultra-rare disorders (PMC11537879)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11537879/); versión Frontiers Med 2024, DOI 10.3389/fmed.2024.1489094 — [enlace](https://www.frontiersin.org/journals/medicine/articles/10.3389/fmed.2024.1489094/full)
- Ejemplo N-of-1 de cribado sobre fibroblastos del paciente: plataforma de imagen de alto contenido con 2.807 compuestos aprobados (NCATS Pharmaceutical Collection) para MPS IIIB — mismo resultado de búsqueda, [PMC11537879](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11537879/) (detalles de autores/DOI del artículo original: no verificados)
- Ejemplo iPSC: en un paciente de 18 años con síndrome tipo Leigh ultra-raro se montó una plataforma basada en iPSC para evaluar un panel de fármacos antes de elegir tratamiento — [hPSCreg, "Development of iPSC-based clinical trial selection platform for patients with ultrarare diseases"](https://hpscreg.eu/browse/publication/2890)
- Ejemplo de mecanismo -> reposicionamiento: en 3 pacientes con iMCD-TAFRO refractarios a bloqueo de IL-6, proteómica sérica identificó activación de mTOR en los brotes y sirolimus indujo remisión en los 3; JCI 2019;129(10):4451-4463 (DOI 10.1172/JCI126091, deducido de la URL) — [JCI / PMC6763220](https://pmc.ncbi.nlm.nih.gov/articles/PMC6763220); ensayo posterior de sirolimus en iMCD: [NCT03933904](https://clinicaltrials.gov/study/NCT03933904)
- Sirolimus en EHE (ultra-rara): dato preclínico de modelo derivado de paciente con 80% de inhibición del volumen tumoral frente a efecto despreciable de doxorrubicina (presentación EMA) — [EMA, Pantziarka](https://www.ema.europa.eu/en/documents/presentation/presentation-repurposing-old-drugs-new-ultra-rare-indication-example-sirolimus-ehe-p-pantziarka_en.pdf)
- Para un ASO individualizado, la prueba in vitro en células del paciente fue la base para lanzar el estudio N-of-1 en <1 año desde el primer contacto (milasen) — [Kim et al., NEJM 2019;381:1644-52, DOI 10.1056/NEJMoa1813279](https://www.ovid.com/journals/nejm/pdf/10.1056/nejmoa1813279~patient-customized-oligonucleotide-therapy-for-a-rare)
- Directrices de consenso para diseño y pruebas de eficacia in vitro de ASO n-of-1 (Aartsma-Rus et al., Nucleic Acid Ther. 2021; DOI 10.1089/nat.2019.0789 deducido de la URL) — [SAGE](https://journals.sagepub.com/doi/pdf/10.1089/nat.2019.0789)
- La FDA, bajo RDEP (anunciado 3 sep 2025), acepta como evidencia confirmatoria: efecto sobre la fisiopatología directa, datos de modelo no clínico, datos farmacodinámicos, casos/acceso expandido, historia natural y controles externos — [Clinical Leader](https://www.clinicalleader.com/doc/new-fda-program-adds-flexibility-for-ultra-rare-disease-submissions-0001); resumen WSGR: [WSGR](https://www.wsgr.com/en/insights/fda-announces-rare-disease-evidence-principles-to-accelerate-treatments-for-ultra-rare-genetic-diseases.html)

### Inferences
- Para MVA/BUB1B, una escalera creíble sería: (i) ensayo en células del paciente o en líneas isogénicas (p. ej. fibroblastos/linfoblastoides) con lectura de segregación cromosómica/aneuploidía/niveles de BubR1; (ii) modelos de ratón BubR1 hipomórfico ya publicados como respaldo (no buscados aquí; lo cubre otro investigador); (iii) uso en paciente solo vía acceso expandido con biomarcadores. Esto es una inferencia, no una norma.
- Los ejemplos de reposicionamiento con mejor trazabilidad (sirolimus) partieron de un biomarcador mecanístico medido en el propio paciente; conviene exigir lo mismo a cualquier candidato.

### Gaps
- No se encontró una fuente regulatoria que fije una "escalera" jerarquizada; es práctica emergente.
- No se verificaron los artículos originales de MPS IIIB y Leigh-like (autores, DOI, resultado clínico).
- Modelos de ratón BubR1 y organoides específicos de MVA: fuera de alcance de esta búsqueda.

## 2. Diseños N-of-1 y de poblaciones pequeñas

### Takeaway
Un N-of-1 es interpretable (no anecdótico) cuando tiene hipótesis y desenlace prefijados, comparación intra-paciente (aleatorización, cegado y/o retirada, múltiples periodos) y análisis que modele tendencia, arrastre y autocorrelación; en la práctica muchos "N-of-1" publicados no alcanzan este estándar.

### Cited Findings
- Revisión sistemática de 74 N-of-1 aleatorizados (2011-2023): solo 13 (17,6%) en enfermedades raras; 77% con cegado (92,3% en raras); washout en 43,2% (61,5% en raras), mediana de 6 periodos de 14 días; análisis heterogéneo (regresión 23%, t-test 24%, bayesiano 14,9%); casi la mitad no definió qué es un respondedor; 66,2% con desenlaces reportados por el paciente — Hawksworth et al., Trials 2024;25:263, DOI 10.1186/s13063-024-08100-1, [PMC11020886](https://pmc.ncbi.nlm.nih.gov/articles/PMC11020886/)
- Los N-of-1 son una forma de diseño de caso único (SCED) y de crossover; pueden considerarse en condiciones raras con pocos participantes — misma revisión, [PMC11020886](https://pmc.ncbi.nlm.nih.gov/articles/PMC11020886/)
- Los modelos bayesianos incorporan tendencia, arrastre (carryover) y autocorrelación; los modelos multinivel permiten combinar varios N-of-1 (N-of-1 agregados) para estimar efecto medio, heterogeneidad y mejorar la inferencia individual — artículo de Harvard Data Science Review (Schmid/Yang; autoría y DOI no verificados, la descarga dio 403) [HDSR](https://hdsr.mitpress.mit.edu/pub/b6efwlql/download/pdf); otras revisiones: [PMC12930012](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12930012/), [PMC10817775](https://pmc.ncbi.nlm.nih.gov/articles/PMC10817775), [PMC8032375](https://pmc.ncbi.nlm.nih.gov/articles/PMC8032375)
- La guía EMA de ensayos en poblaciones pequeñas cubre aleatorización, grupo control, sustitutos, registros, controles históricos/externos y series de casos cuando un ECA no es viable — [EMA guideline (CHMP/EWP/83561/2005)](https://www.ema.europa.eu/en/documents/scientific-guideline/guideline-clinical-trials-small-populations_en.pdf) (el código de referencia es de memoria)
- Guía FDA ASO individualizados (clínica): solo aplica cuando las variantes son únicas (típicamente 1-2 pacientes conocidos) — [FDA, Dec 2021, docket FDA-2021-D-1140](https://www.fda.gov/regulatory-information/search-fda-guidance-documents/ind-submissions-individualized-antisense-oligonucleotide-drug-products-severely-debilitating-or-life)
- Por el criterio "plausible mechanism" de la FDA (feb 2026), la evidencia exige historia natural bien caracterizada, confirmación de acción diana en el paciente y mejora del curso clínico — [Arnold & Porter](https://www.arnoldporter.com/en/perspectives/advisories/2026/02/fda-advances-a-plausible-mechanism-framework-for-rare-disease-drug-development-and-shifts-to)

### Inferences
- Hay dos clases de diseño según la enfermedad: sintomática y reversible (se pueden hacer multiple-baseline/retirada aleatorizada/ABAB) frente a progresiva o de riesgo oncológico, donde la retirada no es ética y el comparador debe ser la historia natural o el propio paciente antes/después con biomarcadores repetidos. MVA con predisposición a cáncer cae en la segunda: por tanto el contrafactual y el biomarcador de aneuploidía serían centrales (inferencia).
- Elementos que distinguen interpretable de anecdótico (síntesis de las fuentes): hipótesis y endpoint prefijados, línea base repetida, cegado cuando sea viable, regla de respondedor definida, análisis bayesiano/serie temporal, y publicación de resultados nulos.

### Gaps
- No obtuve fuentes primarias sobre multiple-baseline y randomised withdrawal (SCED, SCRIBE/CENT 2015) en esta sesión; los nombres son conocimiento general, sin cita.
- Autoría/DOI del artículo HDSR no verificados.

## 3. Rutas regulatorias (FDA y EMA)

### Takeaway
En EE. UU. existe un camino explícito para IND individual (acceso expandido) y desde 2021-2026 guías para terapias individualizadas (ASO 2021; "plausible mechanism" feb 2026; RDEP sep 2025). En la UE no hay guías específicas para N-of-1; se usa uso compasivo (art. 83) o paciente nominal (nacional), más designación huérfana como incentivo.

### Cited Findings
- IND de paciente único / emergencia: solicitud escrita, el tratamiento puede empezar a los 30 días o antes si la FDA lo notifica; en emergencia puede autorizarse por teléfono — [FDA, Expanded Access Part 2 (Sept 2019)](https://www.fda.gov/drugs/information-healthcare-professionals-drugs/expanded-access-part-2-how-submit-single-patient-ind-september-2019)
- Documentos FDA de acceso expandido: [Expanded Access ... Q&A](https://www.fda.gov/regulatory-information/search-fda-guidance-documents/expanded-access-investigational-drugs-treatment-use-questions-and-answers); [IRB Review of Individual Patient EA Submissions](https://www.fda.gov/regulatory-information/search-fda-guidance-documents/institutional-review-board-irb-review-individual-patient-expanded-access-submissions-investigational)
- Serie FDA de borradores sobre ASO individualizados: nonclínico (abr 2021), administrativo/procedimental (ene 2021), clínico y CMC (dic 2021) — [nonclínico](https://www.fda.gov/regulatory-information/search-fda-guidance-documents/nonclinical-testing-individualized-antisense-oligonucleotide-drug-products-severely-debilitating-or), [clínico](https://www.fda.gov/regulatory-information/search-fda-guidance-documents/ind-submissions-individualized-antisense-oligonucleotide-drug-products-severely-debilitating-or-life), [PolicyMed](https://policymed.com/2021/02/fda-releases-draft-guidance-on-ind-submissions-for-individualized-antisense-oligonucleotide-drug-products.html), [RAPS](https://www.raps.org/resource/fda-rolls-out-more-guidance-on-n-of-1-gene-therapi.html)
- Borrador FDA, 23 feb 2026: "Considerations for the Use of Plausible Mechanism Framework to Develop Individualized Therapies that Target Specific Genetic Conditions with Known Biological Cause" (https://www.fda.gov/media/191247/download); cinco elementos (variante identificada, terapia dirigida a la alteración proximal, historia natural, confirmación de target engagement, mejora clínica); evidencia: un estudio adecuado y bien controlado con evidencia confirmatoria; compromisos postcomercialización significativos; se centra en edición génica y terapias de ARN/ASO — [Arnold & Porter](https://www.arnoldporter.com/en/perspectives/advisories/2026/02/fda-advances-a-plausible-mechanism-framework-for-rare-disease-drug-development-and-shifts-to), [AJMC](https://www.ajmc.com/view/fda-proposes-new-approval-pathway-for-individualized-ultra-rare-disease-therapies), [Faegre Drinker](https://www.faegredrinker.com/en/insights/publications/2026/2/fda-publishes-long-awaited-draft-guidance-on-plausible-mechanism-framework-to-spur-the-development-of-new-individualized-medicines)
- Política de un ensayo pivotal por defecto (Makary y Prasad, NEJM, feb 2026; DOI no disponible) — Arnold & Porter, enlace anterior
- RDEP (3 sep 2025): población generalmente <1.000 pacientes en EE. UU., defecto genético de base, progresión rápida, sin alternativas; evidencia sustancial con un estudio adecuado (puede ser de un solo brazo) más evidencia confirmatoria sólida — [Clinical Leader](https://www.clinicalleader.com/doc/new-fda-program-adds-flexibility-for-ultra-rare-disease-submissions-0001), [FDA Law Blog (RDEP y programa "bespoke")](https://www.thefdalawblog.com/2025/09/commissioner-makary-and-center-directors-tidmarsh-and-prasad-announce-2-new-programs-for-ultra-rare-disease-therapies-rdep-and-bespoke). Nota: RDEP está pensado para fármacos con defecto genético, no obviamente para repurposing de fármaco aprobado; no verificado.
- Taller FDA/Duke-Margolis RISE (20 nov 2025) sobre terapias individualizadas — [FDA RISE 2](https://www.fda.gov/industry/rare-disease-innovation-science-and-exploration-rise-workshop-series/rise-workshop-2-individualized-therapies-rise), [informe Duke](https://healthpolicy.duke.edu/sites/default/files/2026-02/Report_Individualized%20Therapies%20on%20the%20RISE.pdf)
- EMA, uso compasivo: guía CHMP adoptada el 19 jul 2007 (art. 83 Reg. 726/2004); para pacientes con enfermedad crónica/gravemente debilitante o amenazante para la vida que no puedan tratarse satisfactoriamente con productos autorizados; distinto del "named patient" (médico obtiene el fármaco del fabricante bajo su responsabilidad) — [EMA guideline](https://www.ema.europa.eu/en/documents/regulatory-procedural-guideline/guideline-compassionate-use-medicinal-products-pursuant-article-83-regulation-ec-no-7262004_en.pdf), [EMA Compassionate use](https://www.ema.europa.eu/en/human-regulatory-overview/research-development/compassionate-use), [Jones Day](https://jonesday.com/de/insights/2010/08/compassionate-use-in-europe-a-patchy-framework-for-early-market-entry) (el PDF de EMA no se pudo leer; la fecha viene del resumen de búsqueda)
- EMA no ha emitido guía específica para N-of-1 (según resumen de búsqueda de una revisión; EMA tiene taller de rutas actuales para terapias individualizadas) — [Toward an Extensible Regulatory Framework for N-of-1 to N-of-Few RNA Therapy, PMC12018499](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12018499/), [EMA presentación Reischl-Kok](https://www.ema.europa.eu/en/documents/presentation/presentation-session-1-current-regulatory-routes-individualised-therapies-i-reischl-kok_en.pdf)
- Orphan: UE prevalencia <=5/10.000 + condición grave + sin tratamiento satisfactorio o beneficio significativo; EE. UU. <200.000 personas — [EMA COMP](https://www.ema.europa.eu/en/documents/presentation/orphan-designation-key-concepts-and-evaluation-criteria-jordi-llinares_en.pdf), [FDA](https://fda.gov/media/109496/download)

### Inferences
- Para un fármaco ya aprobado reposicionado en un único niño (UE vs. EE. UU.; el país del paciente no se conoce ni se asume), la ruta realista es acceso expandido/uso fuera de indicación (off-label) bajo responsabilidad médica; las guías de "plausible mechanism" y RDEP sirven sobre todo como marco de evidencia a imitar, y la designación huérfana solo tendría sentido si hay un patrocinador que busque indicación formal (MVA tiene <50 casos).
- En la UE la vía aplicable para un menor sería la de la autoridad nacional del país correspondiente (uso compasivo/medicamento en situaciones especiales); no buscada aquí.

### Gaps
- No se leyó el texto del borrador de feb 2026 directamente en fda.gov (solo vía resúmenes jurídicos); verificar URL y fecha.
- Equivalentes EMA 2021-2026 concretos (p. ej. posibles papers de reflexión de EMA/ITF para terapias individualizadas) no confirmados; autoridades nacionales no consultadas.
- Estado final (no borrador) de las guías ASO 2021 no verificado.

## 4. Precedentes documentados

### Takeaway
Los precedentes más fuertes combinan mecanismo genético conocido + biomarcador/historia natural. Milasen/atipeksen son terapias nuevas (no reposicionadas); lonafarnib es el caso de fármaco previamente desarrollado (en cáncer) aprobado en ultra-rara sobre un control de historia natural; sirolimus-iMCD es reposicionamiento guiado por proteómica.

### Cited Findings
- Milasen (CLN7/MFSD8): variante SVA retrotransposón + c.1102G>C; ASO modulador de splicing; N-of-1 lanzado en <1 año; sin eventos adversos graves; reducción objetiva de convulsiones por EEG y reporte parental — Kim et al., NEJM 2019;381(17):1644-52, DOI 10.1056/NEJMoa1813279 ([PMID 31597037, de memoria, sin verificar]) — [Harvard compgen](https://compgen.hms.harvard.edu/publications/patient-customized-oligonucleotide-therapy-rare-genetic-disease)
- Atipeksen (A-T, ATM c.7865C>T): marco de ASO individualizado en 235 individuos; 9% "probablemente" y 6% "posiblemente" susceptibles; un niño tratado con buena tolerabilidad durante 3 años — Yu et al., Nature 2023;619:828-836 (DOI 10.1038/s41586-023-06277-0 de memoria, sin verificar) — [Harvard Lee lab](https://compgen.hms.harvard.edu/leelab/publications/framework-individualized-splice-switching-oligonucleotide-therapy); ensayo posterior de A-T: [NCT07215416](https://clinicaltrials.gov/study/NCT07215416). Nota: una fuente de búsqueda atribuyó el trabajo a Nature Medicine 2020; la ficha de Harvard indica Nature 2023.
- Lonafarnib (progeria, HGPS): Gordon et al., JAMA 2018;319(16):1687-1695, DOI 10.1001/jama.2018.3264: cohorte observacional de 63 pacientes tratados de dos ensayos de un solo brazo (ProLon1 n=27, ProLon2 n=36) vs 103 no tratados de un estudio de historia natural; mortalidad 3,7% vs 33,3% a mediana de 2,2 años; HR 0,12 (IC95% 0,01-0,93; p=0,04) — [AJMC](https://www.ajmc.com/view/lonafarnib-monotherapy-extends-survival-in-progeria-), [PRF](https://www.progeriaresearch.org/first-ever-demonstrate-survival-benefit-in-progeria/)
- Aprobación FDA de Zokinvy (lonafarnib) el 20 nov 2020: basada en dos estudios con 62 niños tratados frente a 81 niños no tratados emparejados; vida ~2,5 años más larga de media (prensa) — [Boston Globe/AP](https://www.bostonglobe.com/2020/11/20/business/fda-approves-first-drug-rare-disease-that-causes-premature-aging-children); ensayo del triple: [NCT00916747](https://clinicaltrials.gov/study/NCT00916747); monoterapia fase 2 (NCT00425607, de memoria, sin verificar)
- Base editing personalizado para deficiencia de CPS1 (paciente KJ): diseño y fabricación en ~6 meses, LNP hepático; primera dosis feb 2025 (6-7 meses de edad); a abr 2025 tres dosis sin efectos graves, mayor tolerancia a proteína, mitad de medicación, sin picos de amonio en infecciones — Musunuru et al., NEJM 2025 (DOI/PMID no verificados) — [CHOP](https://chop.edu/news/worlds-first-patient-treated-personalized-crispr-gene-editing-therapy-childrens-hospital), [NBC](https://www.nbcnews.com/health/health-news/pennsylvania-baby-first-receive-personalized-gene-therapy-rcna207100)
- Sirolimus en iMCD-TAFRO refractaria: JCI 2019 (ver sección 1), [NCT03933904](https://clinicaltrials.gov/study/NCT03933904)

### Inferences
- La aprobación de lonafarnib sustentó eficacia con un solo brazo + historia natural, un precedente del tipo de control que podría usarse en MVA (pero MVA carece de cohorte de historia natural de tamaño comparable: <50 casos).
- Los casos de ASO/base editing muestran que la validación se apoya en un biomarcador de la lesión molecular y en seguridad por escalado de dosis, no en eficacia estadística.

### Gaps
- Lonafarnib vs progeria: JAMA 2018 es una comparación frente a controles externos; los límites de sesgo de este diseño (HR con IC amplio) conviene señalarlos.
- Otros precedentes (p. ej. otros ASO de n-Lorem, atalureno, etc.) no investigados; PMIDs de Gordon y Kim no confirmados salvo el DOI.

## 5. Ética/IRB y biomarcadores

### Takeaway
Un uso en paciente único requiere consentimiento (21 CFR 50) y revisión IRB; en EE. UU. la FDA permite exención de revisión por comité completo (21 CFR 56.105) con aprobación de un único miembro del IRB. Los biomarcadores aceptables son farmacodinámicos de acción diana y medidas clínicas objetivas prefijadas; la FDA exige confirmar target engagement.

### Cited Findings
- Acceso expandido = investigación clínica con consentimiento informado según 21 CFR 50; IRB revisa antes del tratamiento (salvo emergencia); con el Form FDA 3926 casilla 10.b el médico solicita exención de reunión plenaria y basta concurrencia de un miembro del IRB — [FDA, Information for IRBs](https://www.fda.gov/news-events/expanded-access/expanded-access-information-institutional-review-boards-irbs); [guía IRB individual](https://www.fda.gov/regulatory-information/search-fda-guidance-documents/institutional-review-board-irb-review-individual-patient-expanded-access-submissions-investigational)
- Guía de IRB para acceso expandido individual (documento FDA media/138496): [FDA](https://fda.gov/media/138496/download)
- El marco de plausible mechanism exige demostrar target engagement con datos mecanísticos/farmacodinámicos y mejora clínica con endpoint definido — [Arnold & Porter](https://www.arnoldporter.com/en/perspectives/advisories/2026/02/fda-advances-a-plausible-mechanism-framework-for-rare-disease-drug-development-and-shifts-to)
- RDEP admite datos farmacodinámicos y de fisiopatología directa como evidencia confirmatoria — [Clinical Leader](https://www.clinicalleader.com/doc/new-fda-program-adds-flexibility-for-ultra-rare-disease-submissions-0001)
- En N-of-1 aleatorizados, 66% usa desenlaces reportados por el paciente y casi la mitad no define respondedor (riesgo de sesgo) — [Hawksworth 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11020886/)
- Milasen: respuesta evaluada con EEG y registro de convulsiones por los padres — NEJM 2019 (enlace sección 4)

### Inferences
- Biomarcadores candidatos para MVA/BUB1B (inferencia, no de fuente regulatoria): nivel de proteína BubR1, tasa de aneuploidía/mosaicismo en linfocitos (cariotipo/FISH/secuenciación de célula única), tasa de errores de segregación en cultivo, marcadores de senescencia (p16, SASP), crecimiento y masa muscular, vigilancia tumoral (imagen/ctDNA). Habría que fijar endpoint primario y regla de respondedor antes de empezar.
- Gestión ética: consentimiento parental + asentimiento, plan de parada por toxicidad, monitorización independiente, y atención al riesgo oncológico (agentes que alteren mitosis pueden empeorar aneuploidía).

### Gaps
- No se encontró guía EMA equivalente específica para comités de ética en uso individual; en la UE sería el comité de ética y la autoridad nacional del país correspondiente (no buscado).
- No hay biomarcadores validados para MVA localizados en esta búsqueda; lo anterior es inferencia.
- No se buscó la guía FDA de biomarcadores (BEST, "Biomarker Qualification") ni la guía de endpoints en enfermedades raras.
