# r_epsilon empirical validation -- summary

Accepted instances analysed: 900 / 1188 rows in results.csv.

## 1. r_val and r_eps distributions
- r_val: {'n': 893, 'mean': 1.0884658454647256, 'median': 1.0, 'std': 0.29949733337689827, 'min': 1.0, 'max': 3.0, 'p25': 1.0, 'p75': 1.0} (UNREACHED fraction 0.0077777777777777776)
- eps=0.01: r_eps {'n': 887, 'mean': 1.0913190529875987, 'median': 1.0, 'std': 0.30348457238375176, 'min': 1.0, 'max': 3.0, 'p25': 1.0, 'p75': 1.0}, UNREACHED fraction 0.014444444444444444, fraction(r_eps > r_val)=0.010078387458006719, mean extra shells (both finite)=1.0
- eps=0.02: r_eps {'n': 881, 'mean': 1.0896708286038592, 'median': 1.0, 'std': 0.3013530181084348, 'min': 1.0, 'max': 3.0, 'p25': 1.0, 'p75': 1.0}, UNREACHED fraction 0.021111111111111112, fraction(r_eps > r_val)=0.0167973124300112, mean extra shells (both finite)=1.0
- eps=0.05: r_eps {'n': 866, 'mean': 1.0935334872979214, 'median': 1.0, 'std': 0.3068086476554212, 'min': 1.0, 'max': 3.0, 'p25': 1.0, 'p75': 1.0}, UNREACHED fraction 0.03777777777777778, fraction(r_eps > r_val)=0.03919372900335946, mean extra shells (both finite)=1.0
- eps=0.1: r_eps {'n': 833, 'mean': 1.1080432172869148, 'median': 1.0, 'std': 0.33302970827173917, 'min': 1.0, 'max': 3.0, 'p25': 1.0, 'p75': 1.0}, UNREACHED fraction 0.07444444444444444, fraction(r_eps > r_val)=0.08846584546472565, mean extra shells (both finite)=1.0526315789473684
- eps=0.25: r_eps {'n': 718, 'mean': 1.1532033426183843, 'median': 1.0, 'std': 0.4007448300916148, 'min': 1.0, 'max': 3.0, 'p25': 1.0, 'p75': 1.0}, UNREACHED fraction 0.20222222222222222, fraction(r_eps > r_val)=0.24524076147816348, mean extra shells (both finite)=1.0909090909090908
- eps=0.5: r_eps {'n': 539, 'mean': 1.25417439703154, 'median': 1.0, 'std': 0.5104450052240099, 'min': 1.0, 'max': 4.0, 'p25': 1.0, 'p75': 1.0}, UNREACHED fraction 0.4011111111111111, fraction(r_eps > r_val)=0.48488241881298993, mean extra shells (both finite)=1.139240506329114
- eps=1.0: r_eps {'n': 269, 'mean': 1.3159851301115242, 'median': 1.0, 'std': 0.5603212829052353, 'min': 1.0, 'max': 3.0, 'p25': 1.0, 'p75': 2.0}, UNREACHED fraction 0.7011111111111111, fraction(r_eps > r_val)=0.7558790593505039, mean extra shells (both finite)=1.1176470588235294

## 2. r_0 == r_val (Remark 3)
- n=900, fraction equal=1.0, exceptions=0

## 3. staircase shape
- counts: {'step': 757, 'ramp': 136, 'flat': 7}
- among non-trivial: ramp=0.1511111111111111, step=0.8411111111111111, flat=0.0077777777777777776

## 4. certificate conservativeness
- realised/beta_up(k_false): {'n': 153, 'mean': 0.5723723206536149, 'median': 0.999999999999998, 'std': 0.48194358514094815, 'min': 0.0, 'max': 1.0000000000000169, 'p25': 7.753822110887211e-16, 'p75': 1.0}
- comparable rows: 884; excluded because beta_up(k_false) was only a lower bound: 16 (of which 3 had a realised error above that lower bound, which is not a violation -- see the note)
- certificate always holds: True (violations: 0); monotonicity violations: 0

## 5. cost
- incremental states: {'n': 900, 'mean': 2.0922222222222224, 'median': 1.0, 'std': 4.007965588635162, 'min': 1.0, 'max': 65.0, 'p25': 1.0, 'p75': 1.0}
- bisection states: {'n': 900, 'mean': 5.137777777777778, 'median': 1.0, 'std': 12.88966247960804, 'min': 1.0, 'max': 172.0, 'p25': 1.0, 'p75': 2.0}
- incremental cheaper (states) fraction: 0.29555555555555557
- greedy_chain_bound tight fraction: 0.9844444444444445
- hard incremental/bisection disagreements: 0
- states saved below r_val (subsample n=107): {'n': 107, 'mean': 1.2429906542056075, 'median': 1.0, 'std': 0.8449753987188822, 'min': 1.0, 'max': 6.0, 'p25': 1.0, 'p75': 1.0}

## 6. comparison with the incumbent mean_abs_bias radius
- subsample n=50
- fraction with any disagreement: 0.84
- mean # thresholds disagreeing (of 7): 1.52
- direction counts: {'new_larger': 23, 'old_larger': 14, 'equal': 8, 'mixed': 5}
- old mean_abs_bias monotone counts: {'violated_1_2': 10, 'violated_2_6': 5, 'violated_6_12': 3, 'violated_5_12': 2, 'violated_23_59': 2, 'violated_13_30': 2, 'violated_9_30': 2, 'violated_12_30': 2, 'violated_3_12': 1, 'violated_72_216': 1, 'violated_48_126': 1, 'violated_10_26': 1, 'violated_18_59': 1, 'violated_7_26': 1, 'violated_60_128': 1, 'violated_2_12': 1, 'violated_9_36': 1, 'violated_103_241': 1, 'violated_12_59': 1, 'violated_30_59': 1, 'violated_44_100': 1, 'violated_47_100': 1, 'monotone': 1, 'violated_39_128': 1, 'violated_45_106': 1, 'violated_1_6': 1, 'violated_7_12': 1, 'violated_20_59': 1, 'violated_13_59': 1, 'violated_29_126': 1}
