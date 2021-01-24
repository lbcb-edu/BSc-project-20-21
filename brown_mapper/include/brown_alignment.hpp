#ifndef BROWN_ALIGN_H
#include <iostream>
#define BROWN_ALIGN_H
    namespace brown {

            enum AlignmentType {GLOBAL, LOCAL, SEMIGLOBAL};

            enum AlignmentDirection {DELETION, INSERTION, MATCH, MISMATCH, NONE};
            
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