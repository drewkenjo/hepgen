paste result_sigT.txt result_sigTL.txt result_sigTT.txt | awk '{printf "%.2f 0.36 %.2f 0.36 %.2f %.2f %.2f %.2f %.2f %.2f %.2f 0\n",$1,$1,$2,$3,$4,$8,$9,$13,$14}' | sed 's/ /\t/g'
