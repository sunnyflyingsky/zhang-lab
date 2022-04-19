awk '$10>200 && $10<=300{print $0}' extension_autodock_final_revised.txt | awk '$8<=-8.3{print $0}' | wc -l
awk '$10>300 && $10<=400{print $0}' extension_autodock_final_revised.txt | awk '$8<=-9{print $0}' | wc -l
awk '$10>400 {print $0}' extension_autodock_final_revised.txt | awk '$8<=-11{print $0}' | wc -l

#
#
#


2835939

101515
220457
245675