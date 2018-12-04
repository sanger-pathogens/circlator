# Change Log

## [Unreleased](https://github.com/sanger-pathogens/circlator/tree/HEAD)

[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v1.5.5...HEAD)

**Fixed bugs:**

- Infinite recursion in clean.py ? [\#123](https://github.com/sanger-pathogens/circlator/issues/123)

**Closed issues:**

- how to use circlator on computer without internet access [\#108](https://github.com/sanger-pathogens/circlator/issues/108)

**Merged pull requests:**

- fix readme and add to codecov [\#141](https://github.com/sanger-pathogens/circlator/pull/141) ([ssjunnebo](https://github.com/ssjunnebo))
- Include Docker  [\#137](https://github.com/sanger-pathogens/circlator/pull/137) ([ssjunnebo](https://github.com/ssjunnebo))
- edit Dockerfile [\#136](https://github.com/sanger-pathogens/circlator/pull/136) ([ssjunnebo](https://github.com/ssjunnebo))
- add Dockerfile [\#135](https://github.com/sanger-pathogens/circlator/pull/135) ([ssjunnebo](https://github.com/ssjunnebo))

## [v1.5.5](https://github.com/sanger-pathogens/circlator/tree/v1.5.5) (2018-01-31)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v1.5.4...v1.5.5)

**Closed issues:**

- nucmer version detected by "circlator progcheck" but not "circlator all" [\#122](https://github.com/sanger-pathogens/circlator/issues/122)
- Canu version is not being detected properly [\#113](https://github.com/sanger-pathogens/circlator/issues/113)

**Merged pull requests:**

- failsafe to prevent infinite recursion [\#126](https://github.com/sanger-pathogens/circlator/pull/126) ([pathdevg](https://github.com/pathdevg))

## [v1.5.4](https://github.com/sanger-pathogens/circlator/tree/v1.5.4) (2018-01-31)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v1.5.3...v1.5.4)

**Closed issues:**

- Problems running progcheck, error with installation? [\#124](https://github.com/sanger-pathogens/circlator/issues/124)

**Merged pull requests:**

- Add tests for external version string checking [\#125](https://github.com/sanger-pathogens/circlator/pull/125) ([andrewjpage](https://github.com/andrewjpage))

## [v1.5.3](https://github.com/sanger-pathogens/circlator/tree/v1.5.3) (2017-11-09)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v1.5.1...v1.5.3)

**Closed issues:**

- Clarification regarding nanopore data [\#119](https://github.com/sanger-pathogens/circlator/issues/119)
- nucmer version isn't detected correctly [\#117](https://github.com/sanger-pathogens/circlator/issues/117)
- Error when running circlator: can't find canu version [\#111](https://github.com/sanger-pathogens/circlator/issues/111)
- progcheck: found nucmer but couldn't get version [\#110](https://github.com/sanger-pathogens/circlator/issues/110)
- Installation problem: command not found after pip3 install circlator [\#109](https://github.com/sanger-pathogens/circlator/issues/109)
- memory usage [\#106](https://github.com/sanger-pathogens/circlator/issues/106)
- Recompile with -fPIC error git / pip3 installs [\#103](https://github.com/sanger-pathogens/circlator/issues/103)
- fixstart with user supplied gene sequences does not work [\#93](https://github.com/sanger-pathogens/circlator/issues/93)

**Merged pull requests:**

- version bump [\#120](https://github.com/sanger-pathogens/circlator/pull/120) ([ssjunnebo](https://github.com/ssjunnebo))
- fixes \#110, \#117 [\#118](https://github.com/sanger-pathogens/circlator/pull/118) ([spock](https://github.com/spock))
- circlator version [\#112](https://github.com/sanger-pathogens/circlator/pull/112) ([ssjunnebo](https://github.com/ssjunnebo))
- update AUTHORS and LICENSE [\#94](https://github.com/sanger-pathogens/circlator/pull/94) ([ssjunnebo](https://github.com/ssjunnebo))

## [v1.5.1](https://github.com/sanger-pathogens/circlator/tree/v1.5.1) (2017-04-21)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v1.5.0...v1.5.1)

**Closed issues:**

- Circlator with short reads? [\#88](https://github.com/sanger-pathogens/circlator/issues/88)

**Merged pull requests:**

- Canu1.5 [\#89](https://github.com/sanger-pathogens/circlator/pull/89) ([martinghunt](https://github.com/martinghunt))

## [v1.5.0](https://github.com/sanger-pathogens/circlator/tree/v1.5.0) (2017-04-03)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v1.4.1...v1.5.0)

**Closed issues:**

- version of SPAdes ? [\#85](https://github.com/sanger-pathogens/circlator/issues/85)
- It looks like circlator keeps a common k-mer at the ends of a contig after changing start position. What for? [\#82](https://github.com/sanger-pathogens/circlator/issues/82)
- Error running circlator with option  --assemble\_not\_only\_assembler [\#80](https://github.com/sanger-pathogens/circlator/issues/80)

**Merged pull requests:**

- Add canu [\#87](https://github.com/sanger-pathogens/circlator/pull/87) ([martinghunt](https://github.com/martinghunt))
- Included option --CanuCorrectedErrorRate to allow user to modify opti… [\#86](https://github.com/sanger-pathogens/circlator/pull/86) ([NicolaDM](https://github.com/NicolaDM))
- Splitting reads for small contigs when using Canu [\#84](https://github.com/sanger-pathogens/circlator/pull/84) ([NicolaDM](https://github.com/NicolaDM))
- Adding the option to use Canu instead of SPAdes [\#83](https://github.com/sanger-pathogens/circlator/pull/83) ([NicolaDM](https://github.com/NicolaDM))

## [v1.4.1](https://github.com/sanger-pathogens/circlator/tree/v1.4.1) (2017-01-24)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v1.4.0...v1.4.1)

**Closed issues:**

- Circularization of BAC contigs [\#79](https://github.com/sanger-pathogens/circlator/issues/79)
- failure to circularise: cannot use this pair because longer match was found [\#75](https://github.com/sanger-pathogens/circlator/issues/75)
- samtools command error [\#66](https://github.com/sanger-pathogens/circlator/issues/66)

**Merged pull requests:**

- Fix assemble not only assembler [\#81](https://github.com/sanger-pathogens/circlator/pull/81) ([martinghunt](https://github.com/martinghunt))

## [v1.4.0](https://github.com/sanger-pathogens/circlator/tree/v1.4.0) (2016-10-17)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v1.3.1...v1.4.0)

**Closed issues:**

- Error running SPAdes.  [\#76](https://github.com/sanger-pathogens/circlator/issues/76)
- Circlator using a gene evicted from Uniprot 1 year ago [\#73](https://github.com/sanger-pathogens/circlator/issues/73)

**Merged pull requests:**

- Fixstart expose mincluster opt [\#77](https://github.com/sanger-pathogens/circlator/pull/77) ([martinghunt](https://github.com/martinghunt))

## [v1.3.1](https://github.com/sanger-pathogens/circlator/tree/v1.3.1) (2016-09-27)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v1.3.0...v1.3.1)

**Closed issues:**

- SPAdes 3.7.1 recommendation [\#72](https://github.com/sanger-pathogens/circlator/issues/72)
- Spades could not be detected [\#68](https://github.com/sanger-pathogens/circlator/issues/68)
- fixstart fails if spaces in sequence names [\#64](https://github.com/sanger-pathogens/circlator/issues/64)

**Merged pull requests:**

- Update dnaa [\#74](https://github.com/sanger-pathogens/circlator/pull/74) ([martinghunt](https://github.com/martinghunt))

## [v1.3.0](https://github.com/sanger-pathogens/circlator/tree/v1.3.0) (2016-08-12)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v1.2.1...v1.3.0)

**Closed issues:**

- Spades version detection error [\#67](https://github.com/sanger-pathogens/circlator/issues/67)
- minimus2 using hard-coded /usr/local/bin/show-coords  [\#65](https://github.com/sanger-pathogens/circlator/issues/65)
- Temporary folders not cleaned up when process killed [\#58](https://github.com/sanger-pathogens/circlator/issues/58)
- Add --verbose to the other subcommands [\#57](https://github.com/sanger-pathogens/circlator/issues/57)
- Spades version detection bit wonky? [\#56](https://github.com/sanger-pathogens/circlator/issues/56)
- circlator fails and still creates the final output file 06.fixstart.fasta [\#55](https://github.com/sanger-pathogens/circlator/issues/55)
- Error running circlator 0.14.0 to reproduce results from published article [\#54](https://github.com/sanger-pathogens/circlator/issues/54)

**Merged pull requests:**

- Add --debug option to progcheck [\#71](https://github.com/sanger-pathogens/circlator/pull/71) ([martinghunt](https://github.com/martinghunt))
- Spades python2 issues [\#70](https://github.com/sanger-pathogens/circlator/pull/70) ([martinghunt](https://github.com/martinghunt))
- Refactor fixstart [\#69](https://github.com/sanger-pathogens/circlator/pull/69) ([martinghunt](https://github.com/martinghunt))

## [v1.2.1](https://github.com/sanger-pathogens/circlator/tree/v1.2.1) (2016-04-19)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v1.2.0...v1.2.1)

**Closed issues:**

- disabling mismatch correction in SPAdes [\#51](https://github.com/sanger-pathogens/circlator/issues/51)
- ImportError: No module named 'circlator' [\#48](https://github.com/sanger-pathogens/circlator/issues/48)

**Merged pull requests:**

- Bug fix running version task [\#63](https://github.com/sanger-pathogens/circlator/pull/63) ([martinghunt](https://github.com/martinghunt))
- Version bump 1.2.1 [\#62](https://github.com/sanger-pathogens/circlator/pull/62) ([martinghunt](https://github.com/martinghunt))
- Summarise circularizations [\#61](https://github.com/sanger-pathogens/circlator/pull/61) ([martinghunt](https://github.com/martinghunt))
- Version reporting improvements [\#60](https://github.com/sanger-pathogens/circlator/pull/60) ([martinghunt](https://github.com/martinghunt))
- add verbose option [\#59](https://github.com/sanger-pathogens/circlator/pull/59) ([martinghunt](https://github.com/martinghunt))

## [v1.2.0](https://github.com/sanger-pathogens/circlator/tree/v1.2.0) (2016-03-09)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v1.1.5...v1.2.0)

**Merged pull requests:**

- Expose spades careful only assemble [\#53](https://github.com/sanger-pathogens/circlator/pull/53) ([martinghunt](https://github.com/martinghunt))

## [v1.1.5](https://github.com/sanger-pathogens/circlator/tree/v1.1.5) (2016-03-09)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v1.1.4...v1.1.5)

**Merged pull requests:**

- Mapping tests [\#50](https://github.com/sanger-pathogens/circlator/pull/50) ([martinghunt](https://github.com/martinghunt))
- support for samtools \> 1.2 [\#49](https://github.com/sanger-pathogens/circlator/pull/49) ([gedankenstuecke](https://github.com/gedankenstuecke))

## [v1.1.4](https://github.com/sanger-pathogens/circlator/tree/v1.1.4) (2016-02-19)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v1.1.3...v1.1.4)

**Closed issues:**

- too many values to unpack error on merge command [\#44](https://github.com/sanger-pathogens/circlator/issues/44)

**Merged pull requests:**

- Nucmer simplify true [\#47](https://github.com/sanger-pathogens/circlator/pull/47) ([martinghunt](https://github.com/martinghunt))
- Make usage clearer [\#46](https://github.com/sanger-pathogens/circlator/pull/46) ([martinghunt](https://github.com/martinghunt))

## [v1.1.3](https://github.com/sanger-pathogens/circlator/tree/v1.1.3) (2016-01-25)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v1.1.2...v1.1.3)

**Closed issues:**

- contig.fastg went missing because of SPAdes 3.6.1 output changes [\#39](https://github.com/sanger-pathogens/circlator/issues/39)

**Merged pull requests:**

- Bug running merge \(issue 44\) [\#45](https://github.com/sanger-pathogens/circlator/pull/45) ([martinghunt](https://github.com/martinghunt))

## [v1.1.2](https://github.com/sanger-pathogens/circlator/tree/v1.1.2) (2015-11-30)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v1.1.1...v1.1.2)

**Merged pull requests:**

- Spades 3.6.2 [\#43](https://github.com/sanger-pathogens/circlator/pull/43) ([martinghunt](https://github.com/martinghunt))

## [v1.1.1](https://github.com/sanger-pathogens/circlator/tree/v1.1.1) (2015-10-30)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v1.1.0...v1.1.1)

**Merged pull requests:**

- Spades 3.6.1 [\#42](https://github.com/sanger-pathogens/circlator/pull/42) ([martinghunt](https://github.com/martinghunt))

## [v1.1.0](https://github.com/sanger-pathogens/circlator/tree/v1.1.0) (2015-10-20)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v1.0.2...v1.1.0)

**Closed issues:**

- Numpy for python3 [\#30](https://github.com/sanger-pathogens/circlator/issues/30)

**Merged pull requests:**

- More logging; better spades assembly; update dependencies [\#41](https://github.com/sanger-pathogens/circlator/pull/41) ([martinghunt](https://github.com/martinghunt))
- Spades try all kmers [\#40](https://github.com/sanger-pathogens/circlator/pull/40) ([martinghunt](https://github.com/martinghunt))

## [v1.0.2](https://github.com/sanger-pathogens/circlator/tree/v1.0.2) (2015-09-03)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v1.0.1...v1.0.2)

**Closed issues:**

- circlator progcheck report PATH of binary? [\#36](https://github.com/sanger-pathogens/circlator/issues/36)

**Merged pull requests:**

- Dnaa filtering [\#38](https://github.com/sanger-pathogens/circlator/pull/38) ([martinghunt](https://github.com/martinghunt))

## [v1.0.1](https://github.com/sanger-pathogens/circlator/tree/v1.0.1) (2015-09-01)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v1.0.0...v1.0.1)

**Merged pull requests:**

- Progcheck report full paths [\#37](https://github.com/sanger-pathogens/circlator/pull/37) ([martinghunt](https://github.com/martinghunt))

## [v1.0.0](https://github.com/sanger-pathogens/circlator/tree/v1.0.0) (2015-08-24)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v0.16.1...v1.0.0)

**Merged pull requests:**

- Generalise get\_dnaa; write versions etc info file [\#35](https://github.com/sanger-pathogens/circlator/pull/35) ([martinghunt](https://github.com/martinghunt))

## [v0.16.1](https://github.com/sanger-pathogens/circlator/tree/v0.16.1) (2015-08-20)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v0.16.0...v0.16.1)

**Merged pull requests:**

- Require fastaq 3.6.1; version bump [\#34](https://github.com/sanger-pathogens/circlator/pull/34) ([martinghunt](https://github.com/martinghunt))
- Add travis [\#32](https://github.com/sanger-pathogens/circlator/pull/32) ([martinghunt](https://github.com/martinghunt))
- Added AUTHORS with mh12 as author [\#29](https://github.com/sanger-pathogens/circlator/pull/29) ([aslett1](https://github.com/aslett1))

## [v0.16.0](https://github.com/sanger-pathogens/circlator/tree/v0.16.0) (2015-08-12)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v0.15.1...v0.16.0)

**Merged pull requests:**

- New task progcheck [\#28](https://github.com/sanger-pathogens/circlator/pull/28) ([martinghunt](https://github.com/martinghunt))
- New task get\_dnaa; update dnaA fasta file [\#27](https://github.com/sanger-pathogens/circlator/pull/27) ([martinghunt](https://github.com/martinghunt))

## [v0.15.1](https://github.com/sanger-pathogens/circlator/tree/v0.15.1) (2015-07-31)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v0.15.0...v0.15.1)

**Merged pull requests:**

- bug fix fixstart --ignore option [\#26](https://github.com/sanger-pathogens/circlator/pull/26) ([martinghunt](https://github.com/martinghunt))

## [v0.15.0](https://github.com/sanger-pathogens/circlator/tree/v0.15.0) (2015-07-29)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v0.14.2...v0.15.0)

**Merged pull requests:**

- sanity check input; contig filtering [\#25](https://github.com/sanger-pathogens/circlator/pull/25) ([martinghunt](https://github.com/martinghunt))
- Remove mention of bacteria [\#24](https://github.com/sanger-pathogens/circlator/pull/24) ([martinghunt](https://github.com/martinghunt))
- help exits without error [\#23](https://github.com/sanger-pathogens/circlator/pull/23) ([andrewjpage](https://github.com/andrewjpage))

## [v0.14.2](https://github.com/sanger-pathogens/circlator/tree/v0.14.2) (2015-07-10)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v0.14.1...v0.14.2)

**Merged pull requests:**

- add homebrew command [\#22](https://github.com/sanger-pathogens/circlator/pull/22) ([andrewjpage](https://github.com/andrewjpage))
- Fix a few typos [\#21](https://github.com/sanger-pathogens/circlator/pull/21) ([martinghunt](https://github.com/martinghunt))

## [v0.14.1](https://github.com/sanger-pathogens/circlator/tree/v0.14.1) (2015-06-24)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v0.14.0...v0.14.1)

**Merged pull requests:**

- do not ignore log files [\#20](https://github.com/sanger-pathogens/circlator/pull/20) ([martinghunt](https://github.com/martinghunt))
- Clean refactor [\#19](https://github.com/sanger-pathogens/circlator/pull/19) ([martinghunt](https://github.com/martinghunt))

## [v0.14.0](https://github.com/sanger-pathogens/circlator/tree/v0.14.0) (2015-06-04)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v0.13.2...v0.14.0)

**Closed issues:**

- Install failed - bio-assembly-refinement/setup.py", line 6, [\#13](https://github.com/sanger-pathogens/circlator/issues/13)

**Merged pull requests:**

- Expose spades kmer [\#18](https://github.com/sanger-pathogens/circlator/pull/18) ([martinghunt](https://github.com/martinghunt))

## [v0.13.2](https://github.com/sanger-pathogens/circlator/tree/v0.13.2) (2015-06-02)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v0.13.1...v0.13.2)

**Merged pull requests:**

- Use second best hits [\#17](https://github.com/sanger-pathogens/circlator/pull/17) ([martinghunt](https://github.com/martinghunt))

## [v0.13.1](https://github.com/sanger-pathogens/circlator/tree/v0.13.1) (2015-05-29)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v0.13.0...v0.13.1)

**Merged pull requests:**

- Minimus2 error typos [\#16](https://github.com/sanger-pathogens/circlator/pull/16) ([martinghunt](https://github.com/martinghunt))

## [v0.13.0](https://github.com/sanger-pathogens/circlator/tree/v0.13.0) (2015-05-28)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v0.12.1...v0.13.0)

**Merged pull requests:**

- Use nucmer diagdiff [\#15](https://github.com/sanger-pathogens/circlator/pull/15) ([martinghunt](https://github.com/martinghunt))
- Add wrapper for minimus2 pipeline [\#14](https://github.com/sanger-pathogens/circlator/pull/14) ([martinghunt](https://github.com/martinghunt))
- Update readme [\#12](https://github.com/sanger-pathogens/circlator/pull/12) ([martinghunt](https://github.com/martinghunt))

## [v0.12.1](https://github.com/sanger-pathogens/circlator/tree/v0.12.1) (2015-05-21)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v0.12.0...v0.12.1)

**Merged pull requests:**

- Fix usage typos [\#11](https://github.com/sanger-pathogens/circlator/pull/11) ([martinghunt](https://github.com/martinghunt))

## [v0.12.0](https://github.com/sanger-pathogens/circlator/tree/v0.12.0) (2015-05-20)
[Full Changelog](https://github.com/sanger-pathogens/circlator/compare/v0.2.0...v0.12.0)

**Merged pull requests:**

-  Write final FASTG file [\#10](https://github.com/sanger-pathogens/circlator/pull/10) ([martinghunt](https://github.com/martinghunt))
- Tidy output and usage [\#9](https://github.com/sanger-pathogens/circlator/pull/9) ([martinghunt](https://github.com/martinghunt))
- add circularised to column to merge log [\#8](https://github.com/sanger-pathogens/circlator/pull/8) ([martinghunt](https://github.com/martinghunt))
- trim reads kept from contig ends [\#7](https://github.com/sanger-pathogens/circlator/pull/7) ([martinghunt](https://github.com/martinghunt))
- Use spades graph [\#6](https://github.com/sanger-pathogens/circlator/pull/6) ([martinghunt](https://github.com/martinghunt))
- Do not break reads [\#5](https://github.com/sanger-pathogens/circlator/pull/5) ([martinghunt](https://github.com/martinghunt))

## [v0.2.0](https://github.com/sanger-pathogens/circlator/tree/v0.2.0) (2015-04-16)
**Merged pull requests:**

- make default merge less permissive [\#4](https://github.com/sanger-pathogens/circlator/pull/4) ([martinghunt](https://github.com/martinghunt))
- Add option to not fix start of some contigs [\#3](https://github.com/sanger-pathogens/circlator/pull/3) ([martinghunt](https://github.com/martinghunt))
- Include default genes fa [\#2](https://github.com/sanger-pathogens/circlator/pull/2) ([martinghunt](https://github.com/martinghunt))
- Initial version [\#1](https://github.com/sanger-pathogens/circlator/pull/1) ([martinghunt](https://github.com/martinghunt))


