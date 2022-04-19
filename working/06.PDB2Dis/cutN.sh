sed -e '/error/d' dock6_results_extend_ss.txt | sed -e '/wrong/d' > dock6_results_extend_ss_clean.txt


awk '{print $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' dock6_results_extend_ss_clean.txt | sort | uniq | awk '$2!~"NNNNNN"{print $0}' | sort | uniq | awk '$3!~"------"{print $0}' | sort | uniq > dock6_results_extend_ss_clean_revised.txt