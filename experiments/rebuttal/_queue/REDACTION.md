# Redaction after registration

The files below are listed in `PREREG_SHA256.txt` with the hash they had when they were frozen.
Before the code was shared, local machine paths, the workstation's host name, the cluster's name and a
personal name were replaced by generic wording (for example each launcher named after the workstation became `launch.sh`, a
local repository path became the repository-relative path plus its commit, and "reported to <name>"
became "reported to the authors"; working labels of the planning notes were dropped from the question lines). No design choice, threshold, statistic or decision rule changed.
`PREREG_SHA256.txt` is kept unchanged as the registration record. For each file whose text changed,
the table gives the registered hash, the hash just before the redaction, and the hash now.
A row whose first two hashes differ was already edited after registration; its own DEVIATIONS.md says why.

| file | registered | before redaction | now |
|---|---|---|---|
| `E10_radius_across_valid_sets/FIDELITY.yaml` | `05fd9c53827b7cb3` | `05fd9c53827b7cb3` | `5c18aa4d1cb74b32` |
| `E11_survival_on_cpdag_only_queries/FIDELITY.yaml` | `eb4a8e8dd6d683a5` | `eb4a8e8dd6d683a5` | `95af2de914239be2` |
| `E12_networks_above_500_nodes/FIDELITY.yaml` | `c659701c35a8953c` | `c659701c35a8953c` | `06106612ca60b9e6` |
| `E13_committed_set_consistency_checks/FIDELITY.yaml` | `15d78ec5299e2b30` | `15d78ec5299e2b30` | `a95e16e3db4f79a8` |
| `E14_data_test_for_invalid_sets/FIDELITY.yaml` | `f284f6de23236e5e` | `f284f6de23236e5e` | `eaf4fcd8ffa10d55` |
| `E15_aid_as_survival_baseline/FIDELITY.yaml` | `cbbb6bb944f31842` | `cbbb6bb944f31842` | `c3ee27454a3a896c` |
| `E16_radius_where_separation_undefined/FIDELITY.yaml` | `8c7e7b9ae2851823` | `8c7e7b9ae2851823` | `b7bf44b6ae6958be` |
| `E17_worked_case_and_validity_by_radius/FIDELITY.yaml` | `17a1b10dc6a8ddd3` | `17a1b10dc6a8ddd3` | `9010a4aa679ea93f` |
| `E18_section7_standardised_models/FIDELITY.yaml` | `bed6c321328f2ac3` | `bed6c321328f2ac3` | `f044ccfed2b9010e` |
| `E19_bias_radius_internal_standardisation/PREREG.md` | `9e988d8ce63c6de3` | `9e988d8ce63c6de3` | `2453314cbfca1a92` |
| `E1_cpdag_check_on_lower_bound_rows/FIDELITY.yaml` | `f875946613405e58` | `f875946613405e58` | `4e4dccd0b74b2747` |
| `E20_uumc_generator/FIDELITY.yaml` | `b346912ca2703fd0` | `b346912ca2703fd0` | `85b6177b40ec4674` |
| `E20_uumc_generator/PREREG.md` | `2678de1386afdb04` | `2678de1386afdb04` | `69971087532ae2be` |
| `E21_walkthrough_toy/FIDELITY.yaml` | `d215b1eb43a70bde` | `d215b1eb43a70bde` | `d0f0687905081c96` |
| `E25_amenability_check/FIDELITY.yaml` | `65ba45805a16590d` | `65ba45805a16590d` | `1e57f23ca06543d9` |
| `E25_amenability_check/PREREG.md` | `0b44556fcaf29df4` | `0b44556fcaf29df4` | `bc53bdb9304a5a89` |
| `E26_taeb_edge_addition_tolerance/FIDELITY.yaml` | `b21df18e79101cc3` | `b21df18e79101cc3` | `d5c3822905847aa5` |
| `E26_taeb_edge_addition_tolerance/PREREG.md` | `78cccfca1cd69395` | `78cccfca1cd69395` | `c8e77abef604348b` |
| `E27_amenability_identity_new_populations/FIDELITY.yaml` | `aedc384815d3b91a` | `aedc384815d3b91a` | `2fd8e7e5fbdb6fbd` |
| `E27_amenability_identity_new_populations/PREREG.md` | `4108c5a1ae5ce335` | `4108c5a1ae5ce335` | `5f238bf92c1b8c3f` |
| `E2_running_time_by_graph_size/FIDELITY.yaml` | `3ce621310cfbef8a` | `3ce621310cfbef8a` | `929a6e7ccb9b6ae1` |
| `E4_class_enumeration_baseline/FIDELITY.yaml` | `864a1e9e309cd0cd` | `864a1e9e309cd0cd` | `feff56b389437dbd` |
| `E4_class_enumeration_baseline/PREREG.md` | `485f4db5050738f9` | `485f4db5050738f9` | `d6a9fdbe58669644` |
| `E5_learned_cpdags_by_sample_size/FIDELITY.yaml` | `c8c7e9172d945dde` | `c8c7e9172d945dde` | `22c4e0d0cbb05f4a` |
| `E6_bootstrap_stability_check/FIDELITY.yaml` | `64a03504995eb47d` | `64a03504995eb47d` | `989050048cdb6407` |
| `E7_one_hidden_variable/FIDELITY.yaml` | `7b0647abfb02e085` | `c07010e5615dc5bb` | `c333952d671872e4` |
| `E8_missed_and_extra_edges/FIDELITY.yaml` | `213f6a7e12086015` | `213f6a7e12086015` | `e8da59f8f414a8a1` |
| `E9_parametric_bootstrap_validity/FIDELITY.yaml` | `03e55fc4d61a840d` | `03e55fc4d61a840d` | `a3a6aca32f8d8aee` |
| `E28_joint_radius_edges/PREREG.md` | `73a8d11bbfd3e47d` | `73a8d11bbfd3e47d` | `1440e15db69672c7` |
| `E29_radius_vs_degree_sep_undefined/PREREG.md` | `c3067b143ff4a57a` | `c3067b143ff4a57a` | `1672f164200795b9` |

Rows for E28 and E29 added on 9 November 2026: the planning-note labels, the workstation's host name and the wording "this cycle / a SUMMARY and a DECISIONS row" were replaced in the two pre-registrations after their runs had finished; the registered hashes in `PREREG_SHA256.txt` are unchanged, and E28's fidelity item `prereg-frozen` checks the run's recorded hash against that file, not against the current text.
