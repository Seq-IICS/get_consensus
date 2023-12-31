#!/bin/bash




while getopts ":p:i:r:l:h" opt; do
  case $opt in
    p)
      arg_path="$OPTARG"
      echo "Option -p with argument: $arg_path"
      ;;
    i)
      arg_interval="$OPTARG"
      echo "Option -i with argument: $arg_interval"
      ;;
    f)
      arg_ref="$OPTARG"
      echo "Option -r with argument: $arg_ref"
      ;;
    l)
      arg_list="$OPTARG"
      echo "Option -l with argument: $arg_list"
      ;;
    \?|h)
      echo "Invalid option: -$OPTARG"
      exit 1
      ;;
  esac
done

cd $arg_path
echo "Exec..."
echo $arg_list
for i in $arg_list; do
  echo "Processing $i"
  minimap2 -a $arg_ref ${i}.gz > ${i}.sam
  samtools view ${i}.sam > ${i}.bam
  samtools sort ${i}.sam > ${i}_sorted.bam
  samtools index ${i}_sorted.bam
  samtools mpileup -d 600000 -A -Q 0 ${i}_sorted.bam | ivar consensus -p ${i}_cns.fasta -q 20  -t 0.5 -n N
done
# for i in (arg_list); do
#   echo "Processing $i"
# done
# for i in $(seq -w 1 "$num"); do
#     echo "Processing $i"
#     minimap2 -a $arg_ref barcode${i}.gz > barcode${i}.sam
#     samtools view barcode${i}.sam > barcode${i}.bam
#     samtools sort barcode${i}.sam > barcode${i}_sorted.bam
#     samtools index barcode${i}_sorted.bam
#     samtools mpileup -d 600000 -A -Q 0 barcode${i}_sorted.bam | ivar consensus -p barcode${i}_cns.fasta -q 20  -t 0.5 -n N
# done