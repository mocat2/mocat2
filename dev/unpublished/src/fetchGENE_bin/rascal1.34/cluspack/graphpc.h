/**************************************************************/
/*                                                            */
/*Procedure de test pour savoir si l'on doit scinder un groupe*/
/*par test d'ajustement de la distribution des similarites    */ 
/*renvoie OUI si on doit separer et NON sinon                 */
/*                                                            */
/**************************************************************/

int test_GRAPHPC(int nb_individus,int nb_dimensions,groupe_t *groupe_pere,groupe_t *groupe1,
		 groupe_t *groupe2);

/**************************************************************************************/
/*                                                                                    */
/*Procedure calculant la probabilite que moins d'un certain pourcentage de transitions*/
/*s'effectuent depuis un groupe de sommets vers un autre groupe de sommets            */
/*                                                                                    */
/**************************************************************************************/

double calcul_probabilite_pourcentage_transitions(int nb_individus,double pb,double pourcentage);
