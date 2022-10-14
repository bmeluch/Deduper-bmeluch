Unit test cases:

UMI not recognized
qname:blah	0	1	5	36	71M	*	0	0	sequence	qualqual

Same coordinates different UMIs
qname:UMI1	0	1	10	36	71M	*	0	0	sequence	qualqual
qname:UMI2	0	1	10	36	71M	*	0	0	sequence	qualqual

Same UMI and position different chromosomes
qname:UMI2	0	1	10	36	71M	*	0	0	sequence	qualqual
qname:UMI2	0	2	10	36	71M	*	0	0	sequence	qualqual

Same UMI and chromosome, different positions
qname:UMI3	0	3	8	36	71M	*	0	0	sequence	qualqual
qname:UMI3	0	3	10	36	71M	*	0	0	sequence	qualqual

Same UMI, chromosome, position, different strands
qname:UMI4	0	4	10	36	71M	*	0	0	sequence	qualqual
qname:UMI4	16	4	10	36	71M	*	0	0	sequence	qualqual

Soft clipped duplicates
qname:UMI8	0	4	38	36	71M	*	0	0	sequence	qualqual
qname:UMI8	0	4	40	36	2S69M	*	0	0	sequence	qualqual

Straightforward duplicates
qname:UMI9	0	4	50	36	71M	*	0	0	sequence	qualqual
qname:UMI9	0	4	50	36	71M	*	0	0	sequence	qualqual