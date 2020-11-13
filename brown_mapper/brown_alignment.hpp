#ifndef BROWN_ALIGN_H
#include <iostream>
#define BROWN_ALIGN_H
    namespace brown {

        enum AlignmentType {GLOBAL, LOCAL, SEMIGLOBAL};

        int bzvz();

        int Align(
            const char* , unsigned int ,
            const char* , unsigned int ,
            AlignmentType ,
            int ,
            int ,
            int ,
            std::string* ,
        unsigned int*);
    }
#endif