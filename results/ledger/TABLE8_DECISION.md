# Table 8 — the decision table, primary supplier (`D_LLM`)

229 certifiable queries, 157 with an exact certificate. The first-failing retraction set is named in the analyst's variables; it is a single claim on 141 of the 157 exact rows. Columns left of the double rule are computed from observables; the audit columns read the truth and say whether the named claim was false, which it was on 68 of 157. First 60 exact rows shown; the full table is `decision_table.csv`.

| network | query | r_claim | holds unless retracted | free rule | r_hop | ‖ named claim false | set at truth |
|---|---|---|---|---|---|---|---|
| Acid_1996 | x1 -> x10 | 1 | x4->x1 | 1 | 1 | yes | invalid |
| Acid_1996 | x1 -> x11 | 1 | x4->x1 | 1 | 1 | yes | invalid |
| Acid_1996 | x1 -> x12 | 1 | x4->x1 | 1 | 1 | yes | invalid |
| Acid_1996 | x1 -> x13 | 1 | x4->x1 | 1 | 1 | yes | invalid |
| Acid_1996 | x1 -> x15 | 1 | x4->x1 | 1 | 1 | yes | invalid |
| Acid_1996 | x1 -> x17 | 1 | x4->x1 | 1 | 1 | yes | invalid |
| Acid_1996 | x1 -> x18 | 1 | x4->x1 | 1 | 1 | yes | invalid |
| Acid_1996 | x1 -> x5 | 1 | x4->x1 | 1 | 1 | yes | invalid |
| Acid_1996 | x1 -> x7 | 1 | x4->x1 | 1 | 1 | yes | invalid |
| Acid_1996 | x1 -> x9 | 1 | x4->x1 | 1 | 1 | yes | invalid |
| Acid_1996 | x4 -> x1 | 1 | x4->x1 | undefined | 1 | yes | invalid |
| Acid_1996 | x4 -> x10 | 1 | x4->x1 | undefined | 1 | yes | invalid |
| Acid_1996 | x4 -> x11 | 1 | x4->x1 | undefined | 1 | yes | invalid |
| Acid_1996 | x4 -> x12 | 1 | x4->x1 | undefined | 1 | yes | invalid |
| Acid_1996 | x4 -> x13 | 1 | x4->x1 | undefined | 1 | yes | invalid |
| Didelez_2010 | Age -> HRT | 1 | Age->Occ | undefined | 1 | no | valid |
| Didelez_2010 | Age -> Occ | 1 | Age->Occ | undefined | 1 | no | valid |
| Didelez_2010 | Age -> S | 1 | Age->Occ | undefined | 1 | no | valid |
| Didelez_2010 | Age -> Smo | 1 | Age->Occ | undefined | 1 | no | valid |
| Didelez_2010 | Age -> TCI | 1 | Age->Occ | undefined | 1 | no | valid |
| Didelez_2010 | Age -> Thist | 1 | Age->Occ | undefined | 1 | no | valid |
| Didelez_2010 | Occ -> HRT | 1 | Age->Occ | 1 | 1 | no | valid |
| Didelez_2010 | Occ -> S | 1 | Age->Occ | 1 | 1 | no | valid |
| Didelez_2010 | Occ -> Smo | 1 | Age->Occ | undefined | 2 | no | valid |
| Didelez_2010 | Occ -> TCI | 1 | Age->Occ | 1 | 1 | no | valid |
| Didelez_2010 | Occ -> Thist | 1 | Age->Occ | undefined | 2 | no | valid |
| Didelez_2010 | Smo -> HRT | 1 | Age->Occ | 2 | 2 | no | valid |
| Didelez_2010 | Smo -> S | 1 | Age->Occ | 2 | 2 | no | valid |
| Didelez_2010 | Smo -> TCI | 1 | Age->Occ | 2 | 2 | no | valid |
| Kampen_2014 | AFF -> AIS | 1 | AFF->ALN | undefined | 1 | no | invalid |
| Kampen_2014 | AFF -> CDR | 1 | AFF->ALN | undefined | 1 | no | invalid |
| Kampen_2014 | AFF -> FTW | 1 | AFF->ALN | undefined | 1 | no | invalid |
| Kampen_2014 | AFF -> SAN | 1 | AFF->ALN | undefined | 1 | no | invalid |
| Kampen_2014 | AIS -> DET | 3 | ALN->SAN and APA->SAN and SAN->AIS | 2 | 3 | yes | invalid |
| Kampen_2014 | AIS -> HOS | 3 | ALN->SAN and APA->SAN and SAN->AIS | 2 | 3 | yes | invalid |
| Kampen_2014 | AIS -> SUS | 3 | ALN->SAN and APA->SAN and SAN->AIS | 2 | 3 | yes | invalid |
| Kampen_2014 | ALN -> PER | 2 | AFF->ALN and ALN->PER | 1 | 1 | no | valid |
| Kampen_2014 | CDR -> DET | 2 | ALN->SAN and APA->SAN | 2 | 3 | yes | valid |
| Polzer_2012 | Caries -> Diabetes | 1 | Caries->Psychosocial | undefined | 1 | yes | invalid |
| Polzer_2012 | Caries -> Mortality | 1 | Caries->Psychosocial | undefined | 1 | yes | invalid |
| Polzer_2012 | Caries -> Psychosocial | 1 | Caries->Psychosocial | undefined | 1 | yes | invalid |
| Polzer_2012 | Caries -> ToothLoss | 1 | Caries->Psychosocial | undefined | 1 | yes | invalid |
| Polzer_2012 | Psychosocial -> Diabetes | 1 | Caries->Psychosocial | undefined | 2 | yes | valid |
| Polzer_2012 | Psychosocial -> Mortality | 1 | Caries->Psychosocial | 1 | 1 | yes | invalid |
| Polzer_2012 | Psychosocial -> Smoking | 1 | Caries->Psychosocial | undefined | 2 | yes | valid |
| Schipf_2010 | A -> PA | 1 | A->PA | undefined | 1 | no | valid |
| Schipf_2010 | A -> S | 1 | A->S | undefined | 1 | no | valid |
| Schipf_2010 | A -> T2DM | 1 | A->PA | undefined | 1 | no | valid |
| Schipf_2010 | A -> TT | 1 | A->PA | undefined | 1 | no | valid |
| Schipf_2010 | A -> U | 1 | A->PA | undefined | 1 | no | valid |
| Schipf_2010 | A -> WC | 1 | A->PA | undefined | 1 | no | valid |
| Schipf_2010 | PA -> T2DM | 1 | A->PA | 1 | 1 | no | valid |
| Schipf_2010 | PA -> TT | 1 | A->PA | 1 | 1 | no | valid |
| Schipf_2010 | PA -> U | 1 | A->PA | 1 | 1 | no | valid |
| Schipf_2010 | PA -> WC | 1 | A->PA | 1 | 1 | no | valid |
| Schipf_2010 | S -> T2DM | 1 | A->S | 2 | 1 | no | valid |
| Schipf_2010 | S -> TT | 1 | A->S | 1 | 1 | no | valid |
| Sebastiani_2005 | ANXA2.5 -> ANXA2.11 | 1 | ANXA2.5->ANXA2.7 | undefined | 1 | yes | invalid |
| Sebastiani_2005 | TGFBR3.8 -> BMP6.10 | 1 | TGFBR3.8->Stroke | undefined | 1 | yes | invalid |
| Shrier_2008 | Coach -> IntraGameProprioception | 1 | Coach->TeamMotivation | undefined | 1 | no | valid |
