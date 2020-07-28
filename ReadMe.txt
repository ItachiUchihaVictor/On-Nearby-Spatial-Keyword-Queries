
Before you use this code, you should know
=======================
	1. This code was used for the empirical study of the TKDE'2019 paper 
	"On Nearby-Fit Spatial Keyword Queries".
	1. This code is developed by Victor Junqiu WEI (wjqjsnj@gmail.com).
	2. This code is written in C/C++.
	3. This code runs under Unix/Linux.
	4. In case that you encounter any problems when using this code,
	please figure out the problem by yourself 
	(The code in fact is easy to read and you can modify it for your own purpose).

Usage
=======================

Step 1: specify the data input information.
<This could be done by editting the "config.txt" file
which format is explained in Appendix I.>

Step 2: compile the source code
make

Step 3: Run the code
./NSKQ

Step 4: Collect the querying results and running statistics 
[you can ignore this step if you don't want to collect the information of
querying results and running statistics]

the querying results are stored in "nskq_result.txt"
which format is explained in Appendix IV.

the running statistics are stored in "nskq_stat.txt"
which format is explained in Appendix V.




Appendix

I. The format of Config.txt
=======================

<Cost indicator>
<Algorithm indicator> 
<# of dimensions>
<# of objects>
<Location file>
<# of keywords>
<Keyword file>
<IR-tree option>
<IR-tree file>
<# of query keywords>
<query set size>
<Percentile lower bound>
<Percentile upper bound>


Explanation of the content in config.txt
-----------------------

<Cost indicator>
	= 1: the MaxSum cost function
	= 2: the Diameter cost function

<Algorithm indicator> 
	= 1: MaxSum-Exact (Dia-Exact)
	= 2: MaxSum-Appro (Dia-Appro)
	= 3: Cao-Exact
	= 4: Cao-Appro1
	= 5: Cao-Appro2
	= 7: NSKQ-Appro1
	= 9: Adapt1-Cao1
	= 10: Adapt1-Cao2
	= 13: Adapt1-Long-Appro
	= 14: Adapt2-Cao1
	= 15: Adapt2-Cao2
	= 16: Adapt2-Long
	= 17: NSKQ-Appro2
	= 18: SKEC+
	= 19: SKEC

<# of dimensions>
	: the number of dimensions of spatial space.

<# of objects>
	: the number of spatial objects.

<Location file>
	: the file containing the locations of the spatial objects,
which format is explained in Appendix II.

<# of keywords>
	: the total number of all possible keywords contained by the objects.

<Keyword file>
	: the file containing the keywords of the spatial objects,
which format is explained in Appendix III.

<IR-tree option>
	= 0: the IR-tree has been built (which will be built and stored in <IR-tree file>)
	= 1: the IR-tree has been built (which is stored in <IR-tree file file>)

<IR-tree file>
	: the file for storing a new (or an existing) IR-tree.

<# of query keywords>
	: the number of keywords in the query (i.e., |q.\psi|).

<query set size>
	: the number of queries that will be performed for the same setting, 
	the average statistics based on which will be used.

<Percentile lower bound>
	: the percentile lower bound that is used for generating queries.

<Percentile upper bound>
	: the percentile upper bound that is used for generating queries.


(See file config.txt in the folder for example)

II. The format of <Location file>
=============================

------------------------
<object ID1>, <1st coordinate>, <2nd coordinate>, ..., <m^th coordinate>
<object ID2>, <1st coordinate>, <2nd coordinate>, ..., <m^th coordinate>
<object ID3>, <1st coordinate>, <2nd coordinate>, ..., <m^th coordinate>
...
<object IDn>, <1st coordinate>, <2nd coordinate>, ..., <m^th coordinate>
------------------------

Note that
	n = # of objects
	m = # of dimensions

(See file running-loc in the folder for example)

III. The format of <Keyword file>
=============================

------------------------
<object ID1>, <1st keyword>, <2nd keyword>, ...
<object ID2>, <1st keyword>, <2nd keyword>, ...
<object ID3>, <1st keyword>, <2nd keyword>, ...
...
<object IDn>, <1st keyword>, <2nd keyword>, ...
------------------------

(See file running-doc in the folder for example)

IV. The format of <nskq_results.txt>
=============================

------------------------
Query #1:
Keywords: <Keywords of the query>
======================
<# of objects in the solution>

<object ID1>: <1st relevant of this object> <2nd relevant of this object> ...
<object ID2>: <1st relevant of this object> <2nd relevant of this object> ...
...
<object IDk>: <1st relevant of this object> <2nd relevant of this object> ...


Query #2:
same format as for Query #1.

...

Query #t:
same format as for Query #1
------------------------

Note that 
	k = the number of objects in the solution
	t = <query set size>

(See file nskq_results.txt in the folder for example)

V. The format of <nskq_stat.txt>
=============================
<Average cost function value>

<the time of building the IR-tree>
<the average time of performing a query>

<Memory usage>
<IR-tree memory usage>
<n_1>
<|P|>
<|O'|>
<|\psi|>

<|O_q|>

<ratio_min>
<ratio_max>
<ratio_average>
<ratio_deviation>


(See file nskq_stat.txt in the folder for example)



