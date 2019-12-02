#!/usr/bin/env bash

run_ffmpeg()
{
  echo "Converting ${1}"
  ffmpeg -i ${1} ${2} ${3} ${4}
}

run_single_series()
{
  FILE_EXTENSION=$1
  TRIAL=$2
  T_START=$3
  Q=$4
  K=$5

  for FILE in $(find . -iname *.${FILE_EXTENSION}); do
    mkdir -p "${HOME}/data/$(dirname ${FILE#./})/$(basename ${FILE} .${FILE_EXTENSION})/${TRIAL}"
    run_ffmpeg "${FILE}" \
      "-ss 00:00:${T_START:-10}.000 -vframes ${Q:-50} -vf select='not(mod(n\,${K:-1}))' -vsync vfr" \
      "${HOME}/data/$(dirname ${FILE#./})/$(basename ${FILE} .${FILE_EXTENSION})/${TRIAL}/img-%04d.png" \
      -hide_banner
  done
}

run_all_series()
{
  run_single_series MOV moon-1 10 50 1
  run_single_series MOV moon-2 10 50 5
  run_single_series MOV moon-3 10 50 10
  run_single_series MOV moon-4 10 50 15
  run_single_series MOV moon-5 10 50 25

  run_single_series AVI sun-1 00 50 1
  run_single_series AVI sun-2 00 50 2
  run_single_series AVI sun-3 00 50 4
  run_single_series AVI sun-4 00 50 8
  run_single_series AVI sun-5 00 50 16
  run_single_series AVI sun-6 00 50 32
}

run_all_series

