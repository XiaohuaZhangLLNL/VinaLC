/* 
 * File:   tokenize.h
 * Author: zhang30
 *
 * Created on February 3, 2012, 2:03 PM
 */

#ifndef VINA_TOKENIZE_H
#define	VINA_TOKENIZE_H

#include <algorithm>
#include <cstdlib>
#include <string>
#include <iterator>
#include <sstream>


template < class ContainerT >
void tokenize(const std::string& str, ContainerT& tokens,
              const std::string& delimiters = " ", const bool trimEmpty = true) {
    typedef ContainerT Base; 
    typedef typename Base::value_type ValueType; 
    typedef typename ValueType::size_type SizeType;
    SizeType pos, lastPos = 0;
    while (true) {
        pos = str.find_first_of(delimiters, lastPos);
        if (pos == std::string::npos) {
            pos = str.length();

            if (pos != lastPos || !trimEmpty)
                tokens.push_back(ValueType(str.data() + lastPos, (SizeType)pos - lastPos));

            break;
        } else {
            if (pos != lastPos || !trimEmpty)
                tokens.push_back(ValueType(str.data() + lastPos, (SizeType)pos - lastPos));
        }

        lastPos = pos + 1;
    }
};


template <typename T1, typename T2>
T1 Sstrm(T2 inValue){
    std::stringstream ss;
    ss << inValue;
    T1 outValue;
    ss >> outValue;
    
    return(outValue);
}

#endif	/* VINA_TOKENIZE_H */

