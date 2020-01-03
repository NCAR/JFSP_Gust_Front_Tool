#!/bin/csh -f
#
#  script to manipulate date strings
#  12/30/19  jfb
#
if ( ${#argv} == 0 ) then
  echo usage: $0 10_digit_date  hour_increment
  echo For example: advance_time.csh 2019080312 24
  exit(1)
endif
#
set sinit = `echo "$1" | cut -c 1-8`
set hh    = `echo "$1" | cut -c 9-10`

date +%Y%m%d%H -d "$sinit $hh + $2 hour"
