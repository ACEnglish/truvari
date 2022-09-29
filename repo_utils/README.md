Truvari Functional Tests
========================

Truvari uses [ssshtests](link) to manage functional tests. To keep everything clean,
each module (or sets of tests) is in its own script under `subtests/`. Then, the 
`truvari_ssshtests.sh` calls each of the subtests. 

To facilitate setup, each subtest should run `source setup_test.sh`, which holds
relevant environment variables.



