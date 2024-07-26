#ifndef TRF_PARAMSET_H
#define TRF_PARAMSET_H

struct trfparamset;

int *paramset_datafile(void);
int paramset_maskedfile(void);
int paramset_flankingsequence(void);
int paramset_HTMLoff(void);
int paramset_redundoff(void);
int paramset_ngs(void);

#endif
