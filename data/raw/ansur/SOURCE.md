# ANSUR II — provenance

* Dataset: 2012 Anthropometric Survey of U.S. Army Personnel (ANSUR II), public release
  (U.S. Army Natick Soldier RD&E Center, 2017). Public domain (U.S. Government work).
* Files: `ANSUR II MALE Public.csv` (4,082 subjects) and `ANSUR II FEMALE Public.csv`
  (1,986 subjects), 108 columns (93 measurements + demographics).
* Official landing pages: https://www.openlab.psu.edu/ansur2/ and
  https://ph.health.mil/topics/workplacehealth/ergo/Pages/Anthropometric-Database.aspx
* Copy used here: the public CSVs as mirrored in
  https://github.com/hkair/anthropometric-stats (data/ansur/), because the official hosts
  were not reachable from the build environment. Row counts (4,082 / 1,986) and the
  108-column layout match the official documentation. Re-download from the official
  page and compare the hashes below before publication-grade use.
* SHA-256
  * MALE   e468eb5787b8ba6c716a86378cf476ba9bd8c57e09493f0fc65b1b8f33e48db8
  * FEMALE 4854f4888146c81690b374d07d044008798ae4a21552afb90823ebd64e04b99f
* Units: millimetres; `weightkg` is in hectograms (0.1 kg).
* Columns used: `waistbreadth`, `waistdepth`, `waistcircumference`, `stature`,
  `weightkg`, `abdominalextensiondepthsitting`, `Gender`.
* Caveat (from the survey documentation): military sample, not representative of
  the civilian or patient population; male and female files are analysed separately.
