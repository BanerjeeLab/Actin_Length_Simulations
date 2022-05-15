clear

mkdir -p data
mkdir -p plot
mkdir -p L_data

rm length.txt

rm data/*
rm plot/*


#-----------------------------------------------------

#****************************************************
#	strating loop for parameters
#****************************************************


# loop for Nc ( capping concentration in nano-molar )
for var0 in 0 0.5 1 5 10 20 50 		 										
do
#		->>>>loop_0<<<<-


# loop for Nfr ( formin concentration in nano-molar )
for var1 in 0 		
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

gfortran growthcode_1.f90			
time ./a.out<para_file.txt 

cat data/pl.txt>>length.txt
cat data/m.txt>>m.txt


#	    ->>>> ensemble loop <<<<-
done


#**************** saving ensemble data in single file  *********************


mv length.txt L_data/L_fr_${var1}_cp_${var0}.txt
mv m.txt L_data/P_fr_${var1}_cp_${var0}.txt


#		->>>>loop_1<<<<-
done


#		->>>>loop_0<<<<-
done


#--------- save code in compressed form ------------


echo 'This script is finished'
echo -e "\a"
echo -e "\a"
echo -e "\a"
echo '---------------End----------------'

