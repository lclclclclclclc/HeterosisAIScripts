touch slim_fitness_debug.out;
for ((i = 0; i < 100; i++)); 
do
if [[ $(( $i % 5 )) -eq 0 ]];
then
	echo $i
fi
slim $1 >> slim_fitness_debug.out; 
done
#touch slim_fitness_debug.results;
#grep ERROR slim_fitness_debug.out | wc -l >> slim_fitness_debug.results
#grep female slim_fitness_debug.out | wc -l >> slim_fitness_debug.results
#grep male slim_fitness_debug.out | wc -l >> slim_fitness_debug.results
#echo -------total/f/m--------- >> slim_fitness_debug.results