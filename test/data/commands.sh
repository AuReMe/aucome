#! /bin/bash

aucome --init=alga -v
### BE CAREFUL Put these 2 GBK files in alga/studied_organims/ each in its own directory.
aucome workflow --run=alga --cpu=4 -v --vv
aucome analysis --run=alga -v
