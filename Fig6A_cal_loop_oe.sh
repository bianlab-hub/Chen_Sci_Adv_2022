## shctrl_1
for i in $(seq 22) X; do
	java -jar juicer_tools_1.19.02.jar dump observed KR shctrl_1.hic $i $i BP 10000 |\
  sed '1d' | sed 's/^/'$i'\t/g' | grep -v NaN >> shctrl_1_observed_KR_10kb.txt
done
awk '{print $1 "\t" $2 "\t" $2+10000 "\t" $1 "\t" $3 "\t" $3+10000 "\t" $4}' shctrl_1_observed_KR_10kb.txt > shctrl_1_observed_KR_10kb.bedpe
bedtools pairtopair -a loops.txt -b shctrl_1_observed_KR_10kb.bedpe | sort -k1,1 -k2,2n | uniq > shctrl_1_loops_observed_KR_10kb.txt

## shctrl_2
for i in $(seq 22) X; do
	java -jar juicer_tools_1.19.02.jar dump observed KR shctrl_2.hic $i $i BP 10000 |\
  sed '1d' | sed 's/^/'$i'\t/g' | grep -v NaN >> shctrl_2_observed_KR_10kb.txt
done
awk '{print $1 "\t" $2 "\t" $2+10000 "\t" $1 "\t" $3 "\t" $3+10000 "\t" $4}' shctrl_2_observed_KR_10kb.txt > shctrl_2_observed_KR_10kb.bedpe
bedtools pairtopair -a loops.txt -b shctrl_2_observed_KR_10kb.bedpe | sort -k1,1 -k2,2n | uniq > shctrl_2_loops_observed_KR_10kb.txt

## shsatb1_1
for i in $(seq 22) X; do
	java -jar juicer_tools_1.19.02.jar dump observed KR shsatb1_1.hic $i $i BP 10000 |\
  sed '1d' | sed 's/^/'$i'\t/g' | grep -v NaN >> shsatb1_1_observed_KR_10kb.txt
done
awk '{print $1 "\t" $2 "\t" $2+10000 "\t" $1 "\t" $3 "\t" $3+10000 "\t" $4}' shsatb1_1_observed_KR_10kb.txt > shsatb1_1_observed_KR_10kb.bedpe
bedtools pairtopair -a loops.txt -b shsatb1_1_observed_KR_10kb.bedpe | sort -k1,1 -k2,2n | uniq > shsatb1_1_loops_observed_KR_10kb.txt

## shsatb1_2
for i in $(seq 22) X; do
	java -jar juicer_tools_1.19.02.jar dump observed KR shsatb1_2.hic $i $i BP 10000 |\
  sed '1d' | sed 's/^/'$i'\t/g' | grep -v NaN >> shsatb1_2_observed_KR_10kb.txt
done
awk '{print $1 "\t" $2 "\t" $2+10000 "\t" $1 "\t" $3 "\t" $3+10000 "\t" $4}' shsatb1_2_observed_KR_10kb.txt > shsatb1_2_observed_KR_10kb.bedpe
bedtools pairtopair -a loops.txt -b shsatb1_2_observed_KR_10kb.bedpe | sort -k1,1 -k2,2n | uniq > shsatb1_2_loops_observed_KR_10kb.txt
