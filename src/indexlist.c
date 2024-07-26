#include <stdlib.h>

#include "indexlist.h"

void FreeList(struct index_list *headptr)
{
    struct index_list *nextptr;
    struct index_list *holdptr;

    nextptr = headptr;
    while (nextptr != NULL) {
        holdptr = nextptr;
        nextptr = nextptr->il_next;
        free(holdptr->il_pattern);
        free(holdptr);
    }

    return;
}
