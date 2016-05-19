  seq1=(char *)malloc(sizeof(char *)*longueur_alignement);
  seq2=(char *)malloc(sizeof(char *)*longueur_alignement);
  for(i=0;i<nb_sequences;i++)
    {
      distances[i][i]=0;
      for(j=i+1;j<nb_sequences;j++)
        {
/* remove gap positions from alignment of 2 sequences */
          for(k=0,length=0;k<longueur_alignement;k++)
            {
            if(sequences[i][k]!='O'||sequences[j][k]!='O')
              {
                 seq1[length]=sequences[i][k];
                 seq2[length]=sequences[j][k];
                 length++;
              }
            }
          distances[i][i]=0;
/* calculate number of identities in a fixed window length along alignment and
count number of times identity is greater than a threshold */
          nb_identites_window=0;
          for(k=0;k<length-window;k++)
            {
               nb_identites=0;
               for(l=k;l<length-window;l++)
                {
                    if(sequences[i][l]==sequences[j][l])
                      {
                        nb_identites++;
                      }
                }
               if(nb_identites>0.3*(float)window) nb_identites_window++;
            }
          if(length==0)
            {
              distances[i][j]=1.0;
            }
          else
            {
              distances[i][j]=1.0-(float)nb_identites_window/(float)length;
            }
          distances[j][i]=distances[i][j];
        }
    }
  free(seq1);free(seq2);
}

