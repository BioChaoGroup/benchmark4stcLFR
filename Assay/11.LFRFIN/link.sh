for i in {1,2,3};do 
	for j in {1,2};do 
		ln -s ../../LFR195M00/split/BC0000000$i/sort.$j.fq LFR195M00_$i/input/rawSeq_$j.fq.gz
	done
done
