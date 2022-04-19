#!/usr/bin/python
import sys
from fuzzywuzzy import fuzz
#conda activate SmartNet
#python region2.py -i 3_all_revised4.txt -j SS_revised_removed.txt > results.txt

for arg in range(len(sys.argv)):
    if sys.argv[arg]=='-i':
        data=sys.argv[arg+1]
    elif sys.argv[arg]=='-j':
        structure=sys.argv[arg+1]




rnapdbee_ss_5O61_A="-((((......((((((((((.....(((..(((((((((((((......((((((.[[.[[[.))).....((((((...(((]]]..]].)))......))))))...)))....(((....)))(((((((......))))))).(((((((((.........))))))))).))...((((((.....(...[..)......))))).).........((.....))...(((((.....((.]...))....)))))....((([...((((((((((........[[[[[[.....))))))))))....(.((((..(.....((((((.]]]]]].)))))).(((((.......((((((.((((.(((...(.((((......)))))[[(...).(((.....]])))...)))...)))))))))))))))......)...))))]....((((((((..[......))))))))....).(((((].[[[[)))))..)))......)).))))))).))...............((....))....))......(..((((....))))............)).....)))))))))).....(.((((..(((((((((..)))).)))))..)))).{.(.(........).).((((((..(.(((.((((((((((..(((...(((....)))...)))...((....))..((((.....))))....(((((((....)))))))....)))))))))))))..((.(...((((((.....((((((((((.((((...((((((......))))))...)[)))...(((.(((((.......A))))))))..).))))))))).(.(((.....)))...)..))))))......).))).((((((....(((((....)))))(((((((((.([.((((((((.....(((....(((.((((....(.((((....))))..)....)))).)))..)).).....))))))))....).]))))))).))..((((((((((......)))))).))))....(((......)))...((((.(((...((([[[...(.((((((.....(.....((..((((((((((.....<.(((((.((([.(((.........)))..))).....(]..((......)).).)))))...)))))))))).))...........]]].)....))))))))))...)))))))((((((((..))))))))....))))))(.((((((((((......((((((((....))).)))))...))))))))))..).....))))))(......((.....((((((.......))))))((((....((((..((..(((((((...((.....)).))))))....[[((((.((.(.(((...(((....)))...)))..)))...(((((........)))))]]((((((((((.(.....((((..((((((((...(..((((.((...((((((...(((((((......)))))))..).)))))...)).)))).)...))))))))....)))).((...(((((...((((..(((..((((..(((((....<<....))))).((..(((....(((((((((((.((.....))))))))).....))))..[.)))...))...(((......))).....))))...)))...))))(..)....)))))..))..).)))))))))))))).))).....((......)).)))).((((.....))))..))))..((((....((((((((.((..(....(.......(((((....(.{..).....)))))((((((...((((....((((..)])))....))))..))))))(...).......(((((((((.].((.........))...(((((((.....(((....))).......))))))).}((((.(((..(((((((..(((...((....((((..((((((...))))))...))))....))...))))))))..)))))(((((.........)))))........((..........(((((.......)))))..)..)....))))...)).))))))).)..))).)))))))).)))))).....))[[..(((((.....}]])))))..((((((((.(((.((......((((((((((..(((((.(((.(((((.((.((((((((((((..(...........(((...((......((.(((((....))))).))...)).............)))..)..)))))))))))).))...(((.((((((......)))))))))....))))).)))..)))))..(((..{{..[)))((((((..........]))))))([((...(((((((.(.((((((.......))))))...).((][.[[[))....)))))))...((((..(((...{..)))..)))).(((....))).))]]].].....(((((((.((..]]]]..))))))))).....)}......))))..))))))......(((((.(.(.((((...((.......))...)))).)...)))).).).....(((((.((((((..(((((((((......))))))..))).(((((.....))))).....)))))..)...))).))....(((((((....))).))))..a.))..)))))))))))[...(((((((...((((..(((((((.....(....)....))))))).(((((...((((...(((((((((((.>>.))))))..)))))...)))).)))))...(((....((((.....>.......))))....))).))))..]....)))))))..(((((...)))))....((((.((((.......)))).....(((.((((.(((....(((((....)))))...)))..))))..)))...))))....)))).."
rnapdbee_seq_5O61_A="UAAGUGUUUAAGGGCGCAUGGUGGAUGCCUUGGCACUGGGAGCCGAUGAAGGACGUAGGAGGCUGCGAUAAGCCUCGGGGAGCUGUCAACCGAGCGUUGAUCCGAGGAUGUCCGAAUGGGGAAACCCGGCACGAGUGAUGUCGUGUCACCAGGCGCUGAAUAUAUAGGCGUCUGGGGGGAACGCGGGGAAGUGAAACAUCUCAGUACCCGUAGGAAGAGAAAACAAAAUGUGAUUCCGUGAGUAGUGGCGAGCGAAAGCGGAGGAUGGCUAAACCGUAUGCAUGUGAUACCGGGUAGGGGUUGUGUGUGCGGGGUUGUGGGACCUAUCUUUCCGGCUCUACCUGGCUGGAGGGCAGUGAGAAAAUGUUGUGGUUAGCGGAAAUGGCUUGGGAUGGCCUGCCGUAGACGGUGAGAGCCCGGUACGUGAAAACCCGACGUCUGUCUUGAUGGUGUUCCCGAGUAGCAGCGGGCCCGUGGAAUCUGCUGUGAAUCUGCCGGGACCACCCGGUAAGCCUGAAUACUUCCCAGUGACCGAUAGCGGAUUAGUACCGUGAGGGAAUGGUGAAAAGUACCCCGGGAGGGGAGUGAAAGAGUACCUGAAACCGUGCGCUUACAAUCCGUCAGAGCCCUCGACGUGUCGUGGGGUGAUGGCGUGCCUUUUGAAGAAUGAGCCUGCGAGUCAGGGACAUGUCGCGAGGUUAACCCGGGUGGGGUAGCCGCAGCGAAAGCGAGUCUGAAUAGGGCGUAUCCACACAAGAGUGUGUGGUGUAGUGGUGUGUUCUGGACCCGAAGCGGAGUGAUCUACCCAUGGCCAGGGUGAAGCGCGGGUAAGACCGCGUGGAGGCCCGAACCCACUUAGGUUGAAGACUGAGGGGAUGAGCUGUGGGUAGGGGUGAAAGGCCAAUCAAACUCCGUGAUAGCUGGUUCUCCCCGAAAUGCAUUUAGGUGCAGCGUCGCAUGUUUCUUGCCGGAGGUAGAGCUACUGGAUGGCCGAUGGGCCCCACAGGGUUACUGACGUCAGCCAAACUCCGAAUGCCGGUAAGUCCAAGAGUGCGGCAGUGAGACGGCGGGGGAUAAGCUCCGUGCGUCGAGAGGGAAACAGCCCAGAUCGCCGGCUAAGGCCCCUAAGCGUGUGCUAAGUGGAAAAGGAUGUGCAGUCGCGAAGACAACCAGGAGGUUGGCUUAGAAGCAGCCACCCUUGAAAGAGUGCGUAAUAGCUCACUGGUCAAGUGAUUGUGCGCCGAUAAUGUAGCGGGGCUCAAGCACACCGCCGAAGCCGCGGCAGCCAACGUGUUGGCUGGGUAGGGGAGCGUCCUGCAUCCGGUGAAGCCGCCGAGUGAUCGAGUGGUGGAGGGUGUGGGAGUGAGAAUGCAGGCAUGAGUAGCGAUUAGGCAAGUGAGAACCUUGCCCGCCGAAAGACCAAGGGUUCCUGGGCCAGGCCAGUCCGCCCAGGGUGAGUCGGGACCUAAGGCGAGGCCGACAGGCGUAGUCGAUGGACAACGGGUUGAUAUUCCCGUACCCGUGUAUGUGCGUCCAUGAUGAAUCAGCGGUACUAACCAUCCAAAACCACCGUGACCGCACCUUUCGGGGUGUGGCGUUGGUGGGGCUGCAUGGGACCUUCGUUGGUAGUAGUCAAGCGAUGGGGUGACGCAGGAAGGUAGCCGUACCGGUCAGUGGUAAUACCGGGGUAAGCCUGUAGGGAGUCAGAUAGGUAAAUCCGUCUGGCAUAUAUCCUGAGAGGUGAUGCAUAGCCGAGUGAGGCGAAUUCGGUGAUCCUAUGCUGCCGAGAAAAGCCUCUAGCGAGGACAUACACGGCCCGUACCCCAAACCAACACAGGUGGUCAGGUAGAGAAUACUAAGGCGUACGAGUGAACUAUGGUUAAGGAACUCGGCAAAAUGCCCCCGUAACUUCGGGAGAAGGGGGACCCACAUGGCGUGUAAGCCUUUACGGCCCAAGCGUGAGUGGGUGGCACAAACCAGUGAGAAGCGACUGUUUACUAAAAACACAGGUCCGUGCGAAGUCGCAAGACGAUGUAUACGGACUGACGCCUGCCCGGUGCUGGAAGGUUAAGAGGACCCGUUAACUCCCUUUGGGGGUGAAGCGGAGAAUUUAAGCCCCAGUAAACGGCGGUGGUAACUAUAACCAUCCUAAGGUAGCGAAAUUCCUUGUCGGGUAAGUUCCGACCUGCACGAAUGGCGUAACGACUUCUCAACUGUCUCAACCAUAGACUCGGCGAAAUUGCACUACGAGUAAAGAUGCUCGUUACGCGCGGCAGGACGAAAAGACCCCGGGACCUUCACUACAACUUGGUAUUGGUGCUCGAUACGGUUUGUGUAGGAUAGGUGGGAGACUGUGAAGCUCACACGCCAGUGUGGGUGGAGUCGUUGUUGAAAUACCACUCUGAUCGUAUUGGGCCUCUAACCUCGGACCGUAUAUCCGGUUCAGGGACAGUGCCUGGUGGGUAGUUUAACUGGGGCGGUUGCCUCCUAAAAUGUAACGGAGGCGCCCAAAGGUUCCCUCAACCUGGACGGCAAUCAGGUGUUGAGUGUAAGUGCACAAGGGAGCUUGACUGCGAGACGGACAUGUCGAGCAGGGACGAAAGUCGGGACUAGUGAUCCGGCACCUCUGAGUGGAAGGGGUGUCGCUCAACGGAUAAAAGGUACCCCGGGGAUAACAGGCUGAUCUUCCCCAAGAGUCCAUAUCGACGGGAUGGUUUGGCACCUCGAUGUCGGCUCGUCGCAUCCUGGGGCUGGAGCAGGUCCCAAGGGUUGGGCUGUUCGCCCAUUAAAGCGGCACGCGAGCUGGGUUUAGAACGUCGUGAGACAGUUCGGUCUCUAUCCGCCGCGCGCGUCAGAAGCUUGAGGAAACCUGUCCCUAGUACGAGAGGACCGGGACGGACGAACCUCUGGUAUACCAGUUGUCCCACCAGGGGCACGGCUGGAUAGCCACGUUCGGACAGGAUAACCGCUGAAAGCAUCUAAGCGGGAAACCUCUUCCAAGACCAGGCUUCUCACCCUCUAGGAGGGAUAAGGCCCCCCGCAGACCACGGGAUUGAUAGACCAGACCUGGAAGCCUAGUAAUAGGUGCAGGGAACUGGCACUAACCGGCCGAAAACUUAC"


#+/-71bp, gap 0bp
rnapdbee_ss_6YDP_AA="-----------------------------------------------------------------------[.[[[[..(((((..]]]]](((((((((..((((((((..(.(((.(((..((....))......(((.....(((((...(((((...)))))...))).)).(((....((..((((......))).)))......)))..(.....)))).(((....))).....)))))))......(..((...(((((....)))).))).).)).))))))...(.((([[[.....(((.....((.]]])).......)))..))).)..))))))))).....((([[...(.((((...((..............(((((((......((((((.((........))..))))))......).....(.((....)))))))))......))...))))....((((((...((...((((.........))))...))))))))......{...(((((((.....]]..))))))))))).(((........((....))......)))...)))))(((((.(((((((...((..(((((..((((((((((.....((........))..........(((((((..(......).((.(((...((((((.(....((..(((((....))).......))..))......((.....))......).).)))...)).)))))....)))))))..)).)))))))).(...((((.....))))...........(....(((((((......)))))))....)..)..))))).....(((((((.........)))))))......))...)))))))))).))..(.(..((.(...(((.((..(((((..(((((...(((.......)))..))))).....)))))..)).)))....).))...)..)..(((((((((....)))))))))}.......-----------------------------------------------------------------------"
rnapdbee_seq_6YDP_AA="GUUAAUGUAGCUUAAAUUAUCAAAGCAAGGCACUGAAAAUGCCUAGAUGAGCCUCACAGCUCCAUAAACACACAGGUUUGGUCCUGGCCUUUCUAUUAAUUCUUAAUAAAAUUACACAUGCAAGUAUCCGCGCCCCGGUGAGAAUGCCCUCCAGAUCUUAAAGAUCAAAAGGAGCAGGUAUCAAGCACACCUAUAACGGUAGCUCAUAACGCCUUGCUCAACCACACCCCCACGGGAAACAGCAGUGAUAAAAAUUAAGCCAUGAACGAAAGUUUGACUAAGUUAUAUUAAUUAGAGUUGGUAAAUCUCGUGCCAGCCACCGCGGUCAUACGAUUAACCCAAAUUAAUAGAUCCACGGCGUAAAGAGUGUUUAAGAAAAAAAAUCACAAUAGAGUUAAAUUAUAACUAAGCUGUAAAAAGCCCUAGUUAAAAUAAAAUAACCCACGAAAGUGACUCUAAUAAUCCUGACACACGAUAGCUAGGACCCAAACUGGGAUUAGAUACCCCACUAUGCCUAGCCCUAAACCCAAAUAGUUACAUAACAAAACUAUUCGCCAGAGUACUACUCGCAACUGCCUAAAACUCAAAGGACUUGGCGGUGCUUCACAUCCACCUAGAGGAGCCUGUUCUAUAAUCGAUAAACCCCGAUAGACCUUACCAACCCUUGCCAAUUCAGCCUAUAUACCGCCAUCUUCAGCAAACCCUAAAAAGGAACAAUAGUAAGCACAAUCAUAGCACAUAAAAACGUUAGGUCAAGGUGUAGCUUAUGGGUUGGAAAGAAAUGGGCUACAUUUUCUACAUAAGAAUAUCCACCACACGAAAGUUUUUAUGAAACUAAAAACCAAAGGAGGAUUUAGCAGUAAAUCAAGAAUAGAGUGCUUGAUUGAAUAAGGCCAUGAAGCACGCACACACCGCCCGUCACCCUCCUCAAGCAUGUAGUAAUAAAAAUAACCUAUAUUCAAUUACACAACCAUGCAAGAAGAGACAAGUCGUAACAAGGUAAGCAUACUGGAAAGUGUGCUUGGAUUACCAAAGCAUAGCUUAAACUAAAGCACCUAGUUUACACCUAGAAGAUCCCACAAUGUAUGGGUACUUUGAACCA"


#+/-1098bp/1000bp, gap 27bp

rnapdbee_ss_6YDP_BA="------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------((.(...(..((((.......................................(...[..)........(..(.....((.]...))....)..)......)).)).............((....))....)..............)....))...........(((.......))).{.(.(........).)..(((((........(((((.(..((....))...).)))))..((.(...((((.(.....(((.((((((.(((.....[)))...(((..((((.......<)))).)))..).))))).))).(.(((.....)))...)..).))))......).))..(((((.....((((....))))...(((((((((((.....((........)).....)))))))))))......(((((((((......)))))))))...((......))....((((.(((...................(.(..(((.........)))..).)....(..(((......)))).......................................)))))))........)))))...(((((......))))).........))))).(......((..((((..))))(((..((((......(.(((((....((....))..)))))......(((..(((((((...........(..)..))))))).......))).)...((......)).)))).....)))..((((....(((((((.((..(....(.......((((((((.].((.........))....((((((.............))))))....(((.(((..(((((....)))..)))))(((((........)))))........((..........(((((.......)))))..)..)....)))....)).))))))..)..))).))))))).)))).)).....)..[..(((((.....}].))))).......(((.(((.((......((((.((((((..((((..(((.((((.(((((.(---------------------------))))))....)))).))).))))..(((..}}..[)))((((((..........]))))))((((((..)((.(.(((...[..))))..)).........))))..........(((....))).....)].....))))..))))))......((((..(.((((((...((.......))...))))))...).)).).)......((((.((((((..((.((((((......))))))...)).(((((.....))))).....)))))..)...))).).....(((((((....))).))))..>.))..)))))).((.(((((((......(....).....)))))))(((((...(.(.....).).))))).(((....(((.................)))....))).)).........((((.......))))...----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
rnapdbee_seq_6YDP_BA="GUUAAUGUAGCUUAAAUUAUCAAAGCAAGGCACUGAAAAUGCCUAGAUGAGCCUCACAGCUCCAUAAACACACAGGUUUGGUCCUGGCCUUUCUAUUAAUUCUUAAUAAAAUUACACAUGCAAGUAUCCGCGCCCCGGUGAGAAUGCCCUCCAGAUCUUAAAGAUCAAAAGGAGCAGGUAUCAAGCACACCUAUAACGGUAGCUCAUAACGCCUUGCUCAACCACACCCCCACGGGAAACAGCAGUGAUAAAAAUUAAGCCAUGAACGAAAGUUUGACUAAGUUAUAUUAAUUAGAGUUGGUAAAUCUCGUGCCAGCCACCGCGGUCAUACGAUUAACCCAAAUUAAUAGAUCCACGGCGUAAAGAGUGUUUAAGAAAAAAAAUCACAAUAGAGUUAAAUUAUAACUAAGCUGUAAAAAGCCCUAGUUAAAAUAAAAUAACCCACGAAAGUGACUCUAAUAAUCCUGACACACGAUAGCUAGGACCCAAACUGGGAUUAGAUACCCCACUAUGCCUAGCCCUAAACCCAAAUAGUUACAUAACAAAACUAUUCGCCAGAGUACUACUCGCAACUGCCUAAAACUCAAAGGACUUGGCGGUGCUUCACAUCCACCUAGAGGAGCCUGUUCUAUAAUCGAUAAACCCCGAUAGACCUUACCAACCCUUGCCAAUUCAGCCUAUAUACCGCCAUCUUCAGCAAACCCUAAAAAGGAACAAUAGUAAGCACAAUCAUAGCACAUAAAAACGUUAGGUCAAGGUGUAGCUUAUGGGUUGGAAAGAAAUGGGCUACAUUUUCUACAUAAGAAUAUCCACCACACGAAAGUUUUUAUGAAACUAAAAACCAAAGGAGGAUUUAGCAGUAAAUCAAGAAUAGAGUGCUUGAUUGAAUAAGGCCAUGAAGCACGCACACACCGCCCGUCACCCUCCUCAAGCAUGUAGUAAUAAAAAUAACCUAUAUUCAAUUACACAACCAUGCAAGAAGAGACAAGUCGUAACAAGGUAAGCAUACUGGAAAGUGUGCUUGGAUUACCAAAGCAUAGCUUAAACUAAAGCACCUAGUUUACACCUAGAAGAUCCCACAAUGUAUGGGUACUUUGAACCAAAGCUAGCUCAACAUACUAAACAAAUACAAAAAUACACCAAAAUAAAAUAAAACAUUCACCUAACAUUAAAGUAUAGGAGAUAGAAAUUUUUAUCCUGACGCUAUAGAGAUAGUACCGUAAGGGAAAGAUGAAAGAAUAAAAUAAAAGUAAAAAAAAGCAAAGAUUACCCCUUCUACCUUUUGCAUAAUGGUUUAACCAGAAAAAAUCUAACAAAGAGAACUUUAGCUAGAUACCCCGAAACCAGACGAGCUACCCAUGAGCAGUUUAAAAGAACCAACUCAUCUAUGUGGCAAAAUAGUGAGAAGACUUGUAGGUAGAGGUGAAAAGCCUAACGAGCCUGGUGAUAGCUGGUUGUCCGAGAAAGAAUUUUAGUUCAACCUUAAAAAUACCCCAAAAACCCUAAAUUCCAAUGUAUUUUUAAGAGAUAGUCUAAAAAGGUACAGCUUUUUAGAAACGGAUACAACCUUGACUAGAGAGUAAAAUCUUAAUACUACCAUAGUAGGCCUAAAAGCAGCCAUCAAUUGAGAAAGCGUUAAAGCUCAACAAAUUCACCAACAUAAUCCCAAAAACUAAUAACAAACUCCUAGCCCAAUACCGGACUAAUCUAUUGAAACAUAGAAGCAAUAAUGUUAAUAUGAGUAACAAGAAGCCUUUCUCCUCGCACACGCUUACAUCAGUAACUAAUAAUAUACUGAUAAUUAACAACCAAUAAACCAAAACAACACUAAAACGUUUAUUAAUUACAUUGUUAACCCAACACAGGAGUGCACCAAGGAAAGAUUAAAAGAAGUAAAAGGAACUCGGCAAACACAAACCCCGCCUGUUUACCAAAAACAUCACCUCUAGCAUUACUAGUAUUAGAGGCAAUGCCUGCCCAGUGACACCAGUUUAACGGCCGCGGUAUUCUGACCGUGCAAAGGUAGCAUAAUCACUUGUUCUCCAAAUAAGGACUUGUAUGAAUGGCCACACGAGGGUUUUACUGUCUCUUACUUCCAAUCAGUGAAAUUAACCUUCCCGUGAAGAGGCGGGAAUAAAAAAAUAAGACGAGAAGACCCUAUGGAGCUUUAAUUAACUAUUCCAAAAGUUAAACAACUCAACCACAAAGGGAUAAAACAUAACUUAACAUGGACUAGCAAUUUCGGUUGGGGUGACCUCGGAGUACAAAAAACCCUCCGAGUGAUUUUAAUCUAGACAAACCAGUCAAAAUAACCAUAACAUCACUUAUUGAUCCAAAAUUUUGAUCAACGGAACAAGUUACCCUAGGGAUAACAGCGCAAUCCUGUUCUAGAGUUCCUAUCGACAAUAGGGUUUACGACCUCGAUGUUGGAUCAGGACACCCAAAUGGUGCAACCGCUAUUAAAGGUUCGUUUGUUCAACGAUUAAAGUCCUACGUGAUCUGAGUUCAGACCGGAGCAAUCCAGGUCGGUUUCUAUCUAUUAUAAAUUUCUCCCAGUACGAAAGGACAAGAGAAAUGGGACCAACCUCACAAACGCGUCUCAGAGAUAAUUAAUGAUUUAAUCUUAACCUAAUUAACUCAUAAUAAAUCCAGCCCUAGAACAGGGCACAUUAGGGUGGCAGAGACCGGUAAUUGCGUAAAACUUAAACCUUUAUUACCAGAGGUUCAACUCCUCUCCCUAAUAACAUGUUCAUAAUUAAUAUUCUAAGCCUAAUUAUUCCUAUCCUACUGGCCGUAGCAUUCCUCACCCUAGUAGAACGAAAAGUGCUAGGUUAUAUGCAACUACGAAAAGGACCCAACGUUGUAGGCCCCUACGGCCUACUCCAACCCAUCGCCGAUGCCCUAAAACUAUUCACCAAAGAACCCCUACGACCAGCCACAUCCUCAAUCUCCAUGUUCAUUAUUGCACCAAUCCUAGCCUUAUCCCUAGCACUAACAAUAUGAGUUCCACUACCAAUACCCUACCCUCUAAUCAACAUAAACCUAGGAGUACUAUUCAUGCUAGCCAUGUCAAGCCUAGCAGUCUACUCUAUCCUAUGAUCAGGAUGAGCAUCCAACUCAAAAUACGCACUCAUCGGGGCCCUACGAGCAGUAGCCCAAACAAUCUCAUAUGAAGUAACACUAGCAAUCAUCCUACUAUCAGUACUCCUAAUAAAUGGAUCAUAUACUCUAUCAACCCUAAUCACAACACAAGAGCACAUUUGAAUAAUCUUUACAUCCUGACCCCUAGCCAUAAUAUGAUUUAUCUCAACCCUAGCAGAAACCAACCGAGCCCCGUUCGACCUUACAGAAGGAGAGUCAGAACUUGUAUCAGGCUUUAACGUAGAAUAUGCAGCCGGACCUUUCGCCAUAUUCUUCAUAGCAGAAUAUGCCAACAUCAUCAUAAUAAAUGCAUUUACAGCAAUUCUCUUCCUAGGAGCAUCCCACGACCCACACACACCAGAACUAUAUACAAUCAACUUCGUACUAAAAACACUCGCAUUAACAAUCACCUUCCUAUGAAUCCGAGCAUCAUACCCACGAUUCCGAUACGACCAACUAAUACAUUUACUAUGAAAAAGCUUCCUGCCCCUAACACUAGCUCUAUGUAUAUGACACAUCUCACUCCCU"


"""
#+/-100bp, gap 27bp
rnapdbee_ss_6YDW_BA="----------------------------------------------------------------------------------------------------((.(...(..((((.......................................(...[..)........(..(.....((.]...))....)..)......)).)).............((....))....)..............)....))...........(((.......))).[.(.(........).)..(((((........(((((.(..((....))...).)))))..((.(...((((.(.....(((.((((((.(((.....[)))...(((..((((.......{)))).)))..).))))).))).(.((.......))...)..).))))......).))..(((((.....((((....))))...(((((((((((.....((........)).....)))))))))))......(((((((((......)))))))))...((......))....((((.(((.....................(..(((.........)))..)......(..(((......)))).......................................)))))))........)))))...(((((......))))).........))))).(......((..((((..))))(((..((((......(.(((((....((....))..)))))......(((..(((((((...........(..)..))))))).......))).)...((......)).)))).....)))..((((....(((((((.((..(....(.......((((((((.].((.........))....((((((.............))))))...((((.(((..(((((....)))..)))))(((((........)))))........((..........(((((.......)))))..)..)....))).)..)).))))).).)..))).))))))).)))).)).....)..{..(((((.....]}.))))).......(((.(((.((......((((.((((((..((((..(((.((((.(((((.(---------------------------))))))....)))).))).))))..(((..]]..[)))((((((..........]))))))((((((..)((.(.(((...[..))))..)).........))))..........(((....))).....)].....))))..))))))......((((..(.((((((...((.......))...))))))...).)).).)......((((.((((((..((.((((((......))))))...)).(((((.....))))).....)))))..)...))).).....(((((((....))).))))..}.))..)))))).((.(((((((......(....).....)))))))(((((...(.(.....).).))))).(((.....((.................)).....))).)).........((((.......))))...----------------------------------------------------------------------------------------------------"
rnapdbee_seq_6YDW_BA="GGUAAGCAUACUGGAAAGUGUGCUUGGAUUACCAAAGCAUAGCUUAAACUAAAGCACCUAGUUUACACCUAGAAGAUCCCACAAUGUAUGGGUACUUUGAACCAAAGCUAGCUCAACAUACUAAACAAAUACAAAAAUACACCAAAAUAAAAUAAAACAUUCACCUAACAUUAAAGUAUAGGAGAUAGAAAUUUUUAUCCUGACGCUAUAGAGAUAGUACCGUAAGGGAAAGAUGAAAGAAUAAAAUAAAAGUAAAAAAAAGCAAAGAUUACCCCUUCUACCUUUUGCAUAAUGGUUUAACCAGAAAAAAUCUAACAAAGAGAACUUUAGCUAGAUACCCCGAAACCAGACGAGCUACCCAUGAGCAGUUUAAAAGAACCAACUCAUCUAUGUGGCAAAAUAGUGAGAAGACUUGUAGGUAGAGGUGAAAAGCCUAACGAGCCUGGUGAUAGCUGGUUGUCCGAGAAAGAAUUUUAGUUCAACCUUAAAAAUACCCCAAAAACCCUAAAUUCCAAUGUAUUUUUAAGAGAUAGUCUAAAAAGGUACAGCUUUUUAGAAACGGAUACAACCUUGACUAGAGAGUAAAAUCUUAAUACUACCAUAGUAGGCCUAAAAGCAGCCAUCAAUUGAGAAAGCGUUAAAGCUCAACAAAUUCACCAACAUAAUCCCAAAAACUAAUAACAAACUCCUAGCCCAAUACCGGACUAAUCUAUUGAAACAUAGAAGCAAUAAUGUUAAUAUGAGUAACAAGAAGCCUUUCUCCUCGCACACGCUUACAUCAGUAACUAAUAAUAUACUGAUAAUUAACAACCAAUAAACCAAAACAACACUAAAACGUUUAUUAAUUACAUUGUUAACCCAACACAGGAGUGCACCAAGGAAAGAUUAAAAGAAGUAAAAGGAACUCGGCAAACACAAACCCCGCCUGUUUACCAAAAACAUCACCUCUAGCAUUACUAGUAUUAGAGGCAAUGCCUGCCCAGUGACACCAGUUUAACGGCCGCGGUAUUCUGACCGUGCAAAGGUAGCAUAAUCACUUGUUCUCCAAAUAAGGACUUGUAUGAAUGGCCACACGAGGGUUUUACUGUCUCUUACUUCCAAUCAGUGAAAUUAACCUUCCCGUGAAGAGGCGGGAAUAAAAAAAUAAGACGAGAAGACCCUAUGGAGCUUUAAUUAACUAUUCCAAAAGUUAAACAACUCAACCACAAAGGGAUAAAACAUAACUUAACAUGGACUAGCAAUUUCGGUUGGGGUGACCUCGGAGUACAAAAAACCCUCCGAGUGAUUUUAAUCUAGACAAACCAGUCAAAAUAACCAUAACAUCACUUAUUGAUCCAAAAUUUUGAUCAACGGAACAAGUUACCCUAGGGAUAACAGCGCAAUCCUGUUCUAGAGUUCCUAUCGACAAUAGGGUUUACGACCUCGAUGUUGGAUCAGGACACCCAAAUGGUGCAACCGCUAUUAAAGGUUCGUUUGUUCAACGAUUAAAGUCCUACGUGAUCUGAGUUCAGACCGGAGCAAUCCAGGUCGGUUUCUAUCUAUUAUAAAUUUCUCCCAGUACGAAAGGACAAGAGAAAUGGGACCAACCUCACAAACGCGUCUCAGAGAUAAUUAAUGAUUUAAUCUUAACCUAAUUAACUCAUAAUAAAUCCAGCCCUAGAACAGGGCACAUUAGGGUGGCAGAGACCGGUAAUUGCGUAAAACUUAAACCUUUAUUACCAGAGGUUCAACUCCUCUCCCUAAUAACAUGUUCAUAAUUAAUAUUCUAAGC"
"""


#+/-1098bp/1000bp, gap 27bp
rnapdbee_ss_6YDW_BA="------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------((.(...(..((((.......................................(...[..)........(..(.....((.]...))....)..)......)).)).............((....))....)..............)....))...........(((.......))).[.(.(........).)..(((((........(((((.(..((....))...).)))))..((.(...((((.(.....(((.((((((.(((.....[)))...(((..((((.......{)))).)))..).))))).))).(.((.......))...)..).))))......).))..(((((.....((((....))))...(((((((((((.....((........)).....)))))))))))......(((((((((......)))))))))...((......))....((((.(((.....................(..(((.........)))..)......(..(((......)))).......................................)))))))........)))))...(((((......))))).........))))).(......((..((((..))))(((..((((......(.(((((....((....))..)))))......(((..(((((((...........(..)..))))))).......))).)...((......)).)))).....)))..((((....(((((((.((..(....(.......((((((((.].((.........))....((((((.............))))))...((((.(((..(((((....)))..)))))(((((........)))))........((..........(((((.......)))))..)..)....))).)..)).))))).).)..))).))))))).)))).)).....)..{..(((((.....]}.))))).......(((.(((.((......((((.((((((..((((..(((.((((.(((((.(---------------------------))))))....)))).))).))))..(((..]]..[)))((((((..........]))))))((((((..)((.(.(((...[..))))..)).........))))..........(((....))).....)].....))))..))))))......((((..(.((((((...((.......))...))))))...).)).).)......((((.((((((..((.((((((......))))))...)).(((((.....))))).....)))))..)...))).).....(((((((....))).))))..}.))..)))))).((.(((((((......(....).....)))))))(((((...(.(.....).).))))).(((.....((.................)).....))).)).........((((.......))))...----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
rnapdbee_seq_6YDW_BA="GUUAAUGUAGCUUAAAUUAUCAAAGCAAGGCACUGAAAAUGCCUAGAUGAGCCUCACAGCUCCAUAAACACACAGGUUUGGUCCUGGCCUUUCUAUUAAUUCUUAAUAAAAUUACACAUGCAAGUAUCCGCGCCCCGGUGAGAAUGCCCUCCAGAUCUUAAAGAUCAAAAGGAGCAGGUAUCAAGCACACCUAUAACGGUAGCUCAUAACGCCUUGCUCAACCACACCCCCACGGGAAACAGCAGUGAUAAAAAUUAAGCCAUGAACGAAAGUUUGACUAAGUUAUAUUAAUUAGAGUUGGUAAAUCUCGUGCCAGCCACCGCGGUCAUACGAUUAACCCAAAUUAAUAGAUCCACGGCGUAAAGAGUGUUUAAGAAAAAAAAUCACAAUAGAGUUAAAUUAUAACUAAGCUGUAAAAAGCCCUAGUUAAAAUAAAAUAACCCACGAAAGUGACUCUAAUAAUCCUGACACACGAUAGCUAGGACCCAAACUGGGAUUAGAUACCCCACUAUGCCUAGCCCUAAACCCAAAUAGUUACAUAACAAAACUAUUCGCCAGAGUACUACUCGCAACUGCCUAAAACUCAAAGGACUUGGCGGUGCUUCACAUCCACCUAGAGGAGCCUGUUCUAUAAUCGAUAAACCCCGAUAGACCUUACCAACCCUUGCCAAUUCAGCCUAUAUACCGCCAUCUUCAGCAAACCCUAAAAAGGAACAAUAGUAAGCACAAUCAUAGCACAUAAAAACGUUAGGUCAAGGUGUAGCUUAUGGGUUGGAAAGAAAUGGGCUACAUUUUCUACAUAAGAAUAUCCACCACACGAAAGUUUUUAUGAAACUAAAAACCAAAGGAGGAUUUAGCAGUAAAUCAAGAAUAGAGUGCUUGAUUGAAUAAGGCCAUGAAGCACGCACACACCGCCCGUCACCCUCCUCAAGCAUGUAGUAAUAAAAAUAACCUAUAUUCAAUUACACAACCAUGCAAGAAGAGACAAGUCGUAACAAGGUAAGCAUACUGGAAAGUGUGCUUGGAUUACCAAAGCAUAGCUUAAACUAAAGCACCUAGUUUACACCUAGAAGAUCCCACAAUGUAUGGGUACUUUGAACCAAAGCUAGCUCAACAUACUAAACAAAUACAAAAAUACACCAAAAUAAAAUAAAACAUUCACCUAACAUUAAAGUAUAGGAGAUAGAAAUUUUUAUCCUGACGCUAUAGAGAUAGUACCGUAAGGGAAAGAUGAAAGAAUAAAAUAAAAGUAAAAAAAAGCAAAGAUUACCCCUUCUACCUUUUGCAUAAUGGUUUAACCAGAAAAAAUCUAACAAAGAGAACUUUAGCUAGAUACCCCGAAACCAGACGAGCUACCCAUGAGCAGUUUAAAAGAACCAACUCAUCUAUGUGGCAAAAUAGUGAGAAGACUUGUAGGUAGAGGUGAAAAGCCUAACGAGCCUGGUGAUAGCUGGUUGUCCGAGAAAGAAUUUUAGUUCAACCUUAAAAAUACCCCAAAAACCCUAAAUUCCAAUGUAUUUUUAAGAGAUAGUCUAAAAAGGUACAGCUUUUUAGAAACGGAUACAACCUUGACUAGAGAGUAAAAUCUUAAUACUACCAUAGUAGGCCUAAAAGCAGCCAUCAAUUGAGAAAGCGUUAAAGCUCAACAAAUUCACCAACAUAAUCCCAAAAACUAAUAACAAACUCCUAGCCCAAUACCGGACUAAUCUAUUGAAACAUAGAAGCAAUAAUGUUAAUAUGAGUAACAAGAAGCCUUUCUCCUCGCACACGCUUACAUCAGUAACUAAUAAUAUACUGAUAAUUAACAACCAAUAAACCAAAACAACACUAAAACGUUUAUUAAUUACAUUGUUAACCCAACACAGGAGUGCACCAAGGAAAGAUUAAAAGAAGUAAAAGGAACUCGGCAAACACAAACCCCGCCUGUUUACCAAAAACAUCACCUCUAGCAUUACUAGUAUUAGAGGCAAUGCCUGCCCAGUGACACCAGUUUAACGGCCGCGGUAUUCUGACCGUGCAAAGGUAGCAUAAUCACUUGUUCUCCAAAUAAGGACUUGUAUGAAUGGCCACACGAGGGUUUUACUGUCUCUUACUUCCAAUCAGUGAAAUUAACCUUCCCGUGAAGAGGCGGGAAUAAAAAAAUAAGACGAGAAGACCCUAUGGAGCUUUAAUUAACUAUUCCAAAAGUUAAACAACUCAACCACAAAGGGAUAAAACAUAACUUAACAUGGACUAGCAAUUUCGGUUGGGGUGACCUCGGAGUACAAAAAACCCUCCGAGUGAUUUUAAUCUAGACAAACCAGUCAAAAUAACCAUAACAUCACUUAUUGAUCCAAAAUUUUGAUCAACGGAACAAGUUACCCUAGGGAUAACAGCGCAAUCCUGUUCUAGAGUUCCUAUCGACAAUAGGGUUUACGACCUCGAUGUUGGAUCAGGACACCCAAAUGGUGCAACCGCUAUUAAAGGUUCGUUUGUUCAACGAUUAAAGUCCUACGUGAUCUGAGUUCAGACCGGAGCAAUCCAGGUCGGUUUCUAUCUAUUAUAAAUUUCUCCCAGUACGAAAGGACAAGAGAAAUGGGACCAACCUCACAAACGCGUCUCAGAGAUAAUUAAUGAUUUAAUCUUAACCUAAUUAACUCAUAAUAAAUCCAGCCCUAGAACAGGGCACAUUAGGGUGGCAGAGACCGGUAAUUGCGUAAAACUUAAACCUUUAUUACCAGAGGUUCAACUCCUCUCCCUAAUAACAUGUUCAUAAUUAAUAUUCUAAGCCUAAUUAUUCCUAUCCUACUGGCCGUAGCAUUCCUCACCCUAGUAGAACGAAAAGUGCUAGGUUAUAUGCAACUACGAAAAGGACCCAACGUUGUAGGCCCCUACGGCCUACUCCAACCCAUCGCCGAUGCCCUAAAACUAUUCACCAAAGAACCCCUACGACCAGCCACAUCCUCAAUCUCCAUGUUCAUUAUUGCACCAAUCCUAGCCUUAUCCCUAGCACUAACAAUAUGAGUUCCACUACCAAUACCCUACCCUCUAAUCAACAUAAACCUAGGAGUACUAUUCAUGCUAGCCAUGUCAAGCCUAGCAGUCUACUCUAUCCUAUGAUCAGGAUGAGCAUCCAACUCAAAAUACGCACUCAUCGGGGCCCUACGAGCAGUAGCCCAAACAAUCUCAUAUGAAGUAACACUAGCAAUCAUCCUACUAUCAGUACUCCUAAUAAAUGGAUCAUAUACUCUAUCAACCCUAAUCACAACACAAGAGCACAUUUGAAUAAUCUUUACAUCCUGACCCCUAGCCAUAAUAUGAUUUAUCUCAACCCUAGCAGAAACCAACCGAGCCCCGUUCGACCUUACAGAAGGAGAGUCAGAACUUGUAUCAGGCUUUAACGUAGAAUAUGCAGCCGGACCUUUCGCCAUAUUCUUCAUAGCAGAAUAUGCCAACAUCAUCAUAAUAAAUGCAUUUACAGCAAUUCUCUUCCUAGGAGCAUCCCACGACCCACACACACCAGAACUAUAUACAAUCAACUUCGUACUAAAAACACUCGCAUUAACAAUCACCUUCCUAUGAAUCCGAGCAUCAUACCCACGAUUCCGAUACGACCAACUAAUACAUUUACUAUGAAAAAGCUUCCUGCCCCUAACACUAGCUCUAUGUAUAUGACACAUCUCACUCCCU"





def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub) # use start += 1 to find overlapping matches
        
        
def fix(full,seq,ss):
    frag_revised=[]
    ss_revised=[]
    gap=0
    pos=0
    loc=0
    seq_frag=seq.upper().split("&")
    ss_frag=ss.split("&")

    if(len(seq)!=len(ss)):
        #print("length error")
        frag_final="error"
        ss_final="error"
    else:
        for i in range(len(seq_frag)):
            #print(full.upper().find(frag[i]))
            #print(seq.upper().find(frag[i]))
            if len(seq_frag[i])>=5:
                pos=full.upper().find(seq_frag[i])
                loc=seq.upper().find(seq_frag[i])
            elif len(seq_frag[i])<5 and len(full) >= 20:
                pos_list=list(find_all(full.upper(),seq_frag[i]))
                pos2 = [j for j in pos_list if j >= pos]
                #print(pos_list,pos,pos2)
                pos=pos2[0]
                loc_list=list(find_all(seq.upper(),seq_frag[i]))
                loc2 = [j for j in loc_list if j >= loc]
                loc=loc2[0]
            if i==0:
                if(pos==loc):
                    #print(frag[i])
                    frag_revised.append(seq_frag[i])
                    ss_revised.append(ss_frag[i])
                elif(pos>loc):
                    num=pos-loc
                    #print("."*num+frag[i])
                    frag_revised.append("."*num+seq_frag[i])
                    ss_revised.append("-"*num+ss_frag[i])
                    gap+=num
                elif (pos<loc):
                    frag_revised.append(seq_frag[i])
                    ss_revised.append(ss_frag[i])
            elif i>0 and i<len(seq_frag):
                if(pos==loc):
                    #print(frag[i])
                    frag_revised.append(seq_frag[i])
                    ss_revised.append(ss_frag[i])
                elif(pos>loc):
                    num=pos-loc+1-gap
                    #print(i,pos,loc,gap,num)
                    #print("."*num+frag[i])
                    frag_revised.append("."*num+seq_frag[i])
                    ss_revised.append("-"*num+ss_frag[i])
                    gap+=num-1
                elif (pos<loc):
                    frag_revised.append(seq_frag[i])
                    ss_revised.append(ss_frag[i])
        frag_final=''.join(frag_revised)
        ss_final=''.join(ss_revised)
        #print(len(ss_final), len(full))
        if(len(ss_final)<len(full)):
            count=len(full)-len(ss_final)
            frag_final=frag_final+"."*count
            ss_final=ss_final+"-"*count
            #print(frag_final,len(frag_final))
            #print(ss_final,len(ss_final))
        else:
            #print(frag_final,len(frag_final))
            #print(ss_final,len(ss_final))
            frag_final=frag_final
            ss_final=ss_final
    return frag_final,ss_final


def read_file(filename):
    ids, names, smiles, seqs, structures, locs, frags, affinitys, name2s, weights = [], [], [], [], [], [], [], [], [], []
    for line in open(filename,'r'):
        id, name, smile, seq, structure, loc, frag, affinity, name2, weight = line.split('\t')
        ids.append(id)
        names.append(name)
        smiles.append(smile)
        seqs.append(seq)
        structures.append(structure)
        locs.append(loc)
        frags.append(frag)
        affinitys.append(affinity)
        name2s.append(name2)
        weights.append(weight.replace('\n',''))

    return ids, names, smiles, seqs, structures, locs, frags, affinitys, name2s, weights



def read_structure(filename):
    ids, strands, seqs, structures= [], [], [], []
    for line in open(filename,'r'):
            id, strand, seq, structure= line.split('\t')
            ids.append(id)
            strands.append(strand)
            seqs.append(seq)
            structures.append(structure.replace('\n',''))
    return ids, strands, seqs, structures





ids, names, smiles, seqs, structures, locs, frags, affinitys, name2s, weights = read_file(data)


s_ids, s_strands, s_seqs, s_structures= read_structure(structure)



#print chromosome
#w=open(Out,'w')
for i in range(0,len(open(data).readlines()),1):
    for j in range(0,len(open(structure).readlines()),1):
        if str(ids[i])=="6N5L" or str(ids[i])=="6N5R": #ID error
            pass
        elif str(ids[i])==str(s_ids[j]).upper() and str(seqs[i]) == str(s_seqs[j]).upper():
        	print(ids[i]+ '\t' + names[i] + '\t' + smiles[i] + '\t' + seqs[i].upper() + '\t' + s_structures[j] + '\t' + locs[i] + '\t' + frags[i] + '\t' + affinitys[i] + '\t' + name2s[i] + '\t' + weights[i])
            #print(ids[i]+ '\t' + nums[i] + '\t' + cats[i] + '\t' + species[i] + '\t' + resolutions[i] + '\t' + dates[i] + '\t' +ligands[i] + '\t' + names[i] + '\t' + weights[i] + '\t' + s_structures[j] + '\t' + seq2s[i].upper() + '\t' + modifieds[i] + '\t' + counts[i] + '\t' + strands[i] + '\t' + locs[i] + '\t' + frags[i] + '\t' + rnas[i] + '\t' + proteins[i] + '\t' + pers[i])
        elif str(ids[i])==str(s_ids[j]).upper() and str(seqs[i]) != str(s_seqs[j]).upper():
            if str(s_seqs[j]).upper().find("&")!=-1:
                if len(seqs[i])-len(s_seqs[j])<=3000:
                    frag_final, ss_final = fix(seqs[i],s_seqs[j],s_structures[j])
                    
                    #print(ids[i]+ '\t' + nums[i] + '\t' + cats[i] + '\t' + species[i] + '\t' + resolutions[i] + '\t' + dates[i] + '\t' + ligands[i] + '\t' + names[i] + '\t' + weights[i] + '\t' + ss_final + '\t' + seq2s[i].upper() + '\t' + modifieds[i] + '\t' + counts[i] + '\t' + strands[i] + '\t' + locs[i] + '\t' + frags[i] + '\t' + rnas[i] + '\t' + proteins[i] + '\t' + pers[i])
                    print(ids[i]+ '\t' + names[i] + '\t' + smiles[i] + '\t' + seqs[i].upper() + '\t' + ss_final + '\t' + locs[i] + '\t' + frags[i] + '\t' + affinitys[i] + '\t' + name2s[i] + '\t' + weights[i])
                #else:
                #    print(ids[i]+ '\t' + strands[i] + '\t' + "too_large")
            else:
                if len(seqs[i])-len(s_seqs[j])<=3000:
                    #if fuzz.ratio(s_seqs[j].upper(), seq2s[i]) >=80:
                    #if str(ids[i])=="5O61" and str(strands[i])=="A":
                    #    print(ids[i]+ '\t' + nums[i] + '\t' + cats[i] + '\t' + species[i] + '\t' + resolutions[i] + '\t' + dates[i] + '\t' + ligands[i] + '\t' + names[i] + '\t' + weights[i] + '\t' + rnapdbee_ss_5O61_A + '\t' + rnapdbee_seq_5O61_A.upper() + '\t' + modifieds[i] + '\t' + counts[i] + '\t' + strands[i] + '\t' + locs[i] + '\t' + frags[i] + '\t' + rnas[i] + '\t' + proteins[i] + '\t' + pers[i])
                    if fuzz.partial_ratio(s_seqs[j].upper(), seqs[i]) >=80: #or substirng
                        #print(ids[i]+ '\t' + nums[i] + '\t' + cats[i] + '\t' + species[i] + '\t' + resolutions[i] + '\t' + dates[i] + '\t' + ligands[i] + '\t' + names[i] + '\t' + weights[i] + '\t' + s_structures[j] + '\t' + s_seqs[j].upper() + '\t' + modifieds[i] + '\t' + counts[i] + '\t' + strands[i] + '\t' + locs[i] + '\t' + frags[i] + '\t' + rnas[i] + '\t' + proteins[i] + '\t' + pers[i])
                        print(ids[i]+ '\t' + names[i] + '\t' + smiles[i] + '\t' + s_seqs[j].upper() + '\t' + s_structures[j] + '\t' + locs[i] + '\t' + frags[i] + '\t' + affinitys[i] + '\t' + name2s[i] + '\t' + weights[i])
                else:
                    #print(ids[i]+ '\t' + nums[i] + '\t' + cats[i] + '\t' + species[i] + '\t' + resolutions[i] + '\t' + dates[i] + '\t' + ligands[i] + '\t' + names[i] + '\t' + weights[i] + '\t' + s_structures[j] + '\t' + s_seqs[j].upper() + '\t' + modifieds[i] + '\t' + counts[i] + '\t' + strands[i] + '\t' + locs[i] + '\t' + frags[i] + '\t' + rnas[i] + '\t' + proteins[i] + '\t' + pers[i])
                    print(ids[i]+ '\t' + names[i] + '\t' + smiles[i] + '\t' + s_seqs[j].upper() + '\t' + s_structures[j] + '\t' + locs[i] + '\t' + frags[i] + '\t' + affinitys[i] + '\t' + name2s[i] + '\t' + weights[i])
                #else:
                #    print(ids[i]+ '\t' + "unknown" + '\t'+ s_seqs[j].upper() + '\t' + seq2s[i]) #5O61

        #elif str(ids[i])==str(s_ids[j]).upper() and str(strands[i]) != str(s_strands[j]).upper() and str(seq2s[i]) != str(s_seqs[j]).upper():
        #    print(ids[i]+ '\t' + "unknown")
        #else:
        #    print("error_unknown")