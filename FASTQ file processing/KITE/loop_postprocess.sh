#!/bin/bash

set -e

folders=(
    '/KITE/HTO/Set-1/2022-03-16-1-1'
    '/KITE/HTO/Set-1/2022-03-16-1-2'
    '/KITE/HTO/Set-1/2022-03-16-1-3'
    '/KITE/HTO/Set-1/2022-03-16-1-4'
    
    '/KITE/HTO/Set-2/2022-03-16-1-5'
    '/KITE/HTO/Set-2/2022-03-16-1-6'
    '/KITE/HTO/Set-2/2022-03-16-1-7'
    '/KITE/HTO/Set-2/2022-03-16-1-8'

    '/KITE/HTO/Set-3/2022-03-16-2-1'
    '/KITE/HTO/Set-3/2022-03-16-2-2'
    '/KITE/HTO/Set-3/2022-03-16-2-3'


    '/KITE/ADT/2022-03-16-1-1'
    '/KITE/ADT/2022-03-16-1-2'
    '/KITE/ADT/2022-03-16-1-3'
    '/KITE/ADT/2022-03-16-1-4'
    
    '/KITE/ADT/2022-03-16-1-5'
    '/KITE/ADT/2022-03-16-1-6'
    '/KITE/ADT/2022-03-16-1-7'
    '/KITE/ADT/2022-03-16-1-8'

    '/KITE/ADT/2022-03-16-2-1'
    '/KITE/ADT/2022-03-16-2-2'
    '/KITE/ADT/2022-03-16-2-3'
)

for folder in ${folders[@]}; do

    echo "Starting $folder"
    postprocess_kite_output.sh $folder

done