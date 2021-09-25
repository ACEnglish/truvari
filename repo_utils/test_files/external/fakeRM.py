#!/bin/env python3
import sys
ans = """  bit   perc perc perc  query     position in query    matching   repeat          position in repeat
score   div. del. ins.  sequence  begin end   (left)   repeat     class/family  begin  end    (left) 

  129   16.6  0.0  5.4  1028          2   176    (0) C ALUJO      SINE/Alu        (53)    259     94  
  122   18.3  0.0  5.4  1029          2   176    (0) C ALUJR      SINE/Alu        (53)    259     94  
  134    9.5  4.1  0.0  1474          6   225    (2) + (CCCCATC)N Simple_repeat      1    229    (0)  
  135    9.9  3.6  0.0  1475          6   226    (2) + (CCCCATC)N Simple_repeat      1    229    (0)  
   64   21.1  0.0  0.0  3352          1   104    (0) C L2A        LINE/L2        (279)   3147   3044  
   97    8.4  0.7  0.0  909           1   139    (1) + (GGGA)N    Simple_repeat      1    140    (0)  
  300   10.4  0.5  0.0  910           1   438    (1) + (AGGG)N    Simple_repeat      1    440    (0)  
  360    9.6  0.4  0.0  912           1   514    (1) + (AGGG)N    Simple_repeat      1    516    (0)  
  117    9.7  1.1  0.6  913           1   179    (1) + (GGGA)N    Simple_repeat      1    180    (0)  """


if __name__ == '__main__':
    fasta_name = sys.argv[-1]
    print('fout would be', fasta_name + ".out")
    with open(fasta_name + ".out", 'w') as fout:
        fout.write(ans)
