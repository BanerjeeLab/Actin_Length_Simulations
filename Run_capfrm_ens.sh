clear

mkdir -p data
mkdir -p plot
mkdir -p L_data

rm length.txt
rm allxx.txt

rm data/*
rm plot/*


#-----------------------------------------------------

#****************************************************
#	strating loop for parameters
#****************************************************

## one micro-molar concentration ~ 600 monomers in 1 micro-m^3
## one nano-molar concentration ~ 6 monomers in 10 micro-m^3

# loop for Nc ( capping concentration in nano-molar )
for var0 in 5 		 										
do
#		->>>>loop_0<<<<-


# loop for Nfr ( formin concentration in nano-molar )
for var1 in 5 		
do
#		->>>>loop_1<<<<-

#------------------- Ensemble loop--------------------------

for ven in {1..10}
do

ran=$RANDOM

printf "%f" $ran>'ran.txt'

rm a.out

printf "%f %f" $var0 $var1>'para_file.txt'

echo 'Nc =' ${var0}
echo 'Nfr =' ${var1}
echo 'running the code'	

gfortran PFCPint_5uM.f90			
time ./a.out<para_file.txt 

cat data/l.txt>>length.txt
cat data/m.txt>>m.txt
cat xx.txt>>allxx.txt

#**************** plot *****************

#gnuplot plot_frm.gp


#**************** saving data  *********************

data_flag=0
if [ "$data_flag" -eq 1 ]; then
cd data
mv l.txt L_fr_${var1}.txt
mv m.txt m_fr_${var1}.txt
cd ..
fi





#	    ->>>> ensemble loop <<<<-
done



#***************************************************

#gnuplot plot_nuc.gp

#**************** saving plot  *********************

plot_flag=0
if [ "$plot_flag" -eq 1 ]; then
cd plot
mv filaments.png F_fr_${var1}_cp_${var0}.png
mv mean_size.png Lt_fr_${var1}_cp_${var0}.png
mv Ldist.png Pn_fr_${var1}_cp_${var0}.png
cd ..
fi



mv length.txt L_data/L_fr_${var1}_cp_${var0}.txt
mv m.txt L_data/P_fr_${var1}_cp_${var0}.txt
#mv state.txt S_data/S_fr_${var1}_cp_${var0}.txt


#		->>>>loop_1<<<<-
done


#		->>>>loop_0<<<<-
done


#--------- save code in compressed form ------------

#tar -czvf compCode.tar.gz code make_codeBulk.sh plot*.gp
#mv compCode.tar.gz result/compCode.tar.gz



echo 'This script is finished'
echo -e "\a"
echo -e "\a"
echo -e "\a"
echo '---------------End----------------'

