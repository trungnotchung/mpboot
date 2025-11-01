#include "mutation.h"
#include <fcntl.h>
#include <iomanip>
#include <iostream>
#include <random>
#include <algorithm>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <signal.h>
// Uses one-hot encoding if base is unambiguous
// A:1,C:2,G:4,T:8

// Convert nuc_id back to IUPAC base
char get_nuc(int8_t nuc_id) {
    char ret = 'N';
    //assert ((nuc_id >= 1) && (nuc_id <= 15));
    switch (nuc_id) {
    case 1:
        ret = 'A';
        break;
    case 2:
        ret = 'C';
        break;
    case 3:
        ret = 'M';
        break;
    case 4:
        ret = 'G';
        break;
    case 5:
        ret = 'R';
        break;
    case 6:
        ret = 'S';
        break;
    case 7:
        ret = 'V';
        break;
    case 8:
        ret = 'T';
        break;
    case 9:
        ret = 'W';
        break;
    case 10:
        ret = 'Y';
        break;
    case 11:
        ret = 'H';
        break;
    case 12:
        ret = 'K';
        break;
    case 13:
        ret = 'D';
        break;
    case 14:
        ret = 'B';
        break;
    default:
        ret = 'N';
        break;
    }
    return ret;
}