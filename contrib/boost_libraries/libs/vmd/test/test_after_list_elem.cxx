
//  (C) Copyright Edward Diener 2011-2015
//  Use, modification and distribution are subject to the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt).

#include <boost/vmd/elem.hpp>
#include <boost/vmd/is_empty.hpp>
#include <boost/detail/lightweight_test.hpp>
#include <boost/preprocessor/seq/elem.hpp>

int main()
  {
  
#if BOOST_PP_VARIADICS

  #define BOOST_VMD_REGISTER_ggh (ggh)
  
  #define ANIDENTIFIER ggh
  #define ANUMBER 249
  #define ASEQ (25)(26)(27)
  #define ATUPLE (0,1,2,3,((a,b))((c,d))((e))((f,g,h)))
  #define ALIST (0,(1,(2,(3,BOOST_PP_NIL))))
  #define ALIST2 (10,(11,(12,(13,BOOST_PP_NIL))))
  #define ALIST3 (100,(101,(102,(103,BOOST_PP_NIL))))
  #define ALIST5 (200,(201,(202,(203,BOOST_PP_NIL))))
  #define ANARRAY (3,(a,b,38))
  #define ANARRAY2 (5,(c,d,133,22,15))
  #define ASEQUENCE ANUMBER ALIST ATUPLE ANIDENTIFIER ANARRAY ASEQ
  #define ASEQUENCE2 ANIDENTIFIER ANARRAY2 ASEQ ALIST2 ATUPLE
  #define ASEQUENCE3 ASEQ ANIDENTIFIER ATUPLE ALIST3
  #define ASEQUENCE4
  #define ASEQUENCE5 ALIST5 ASEQ ATUPLE ANIDENTIFIER

  #define A_LIST_PLUS (mmf,(34,(^^,(!,BOOST_PP_NIL)))) 456
  #define PLUS_ALIST yyt (j,(ii%,BOOST_PP_NIL))
  #define JDATA ggh
  #define KDATA (a,(b,BOOST_PP_NIL)) name
  #define ELISTP BOOST_PP_NIL 231
  #define A_SEQ ((25,BOOST_PP_NIL) 3)((26,BOOST_PP_NIL) 4)((27,BOOST_PP_NIL) 5)
  
  BOOST_TEST(!BOOST_VMD_IS_EMPTY(BOOST_VMD_ELEM(1,ASEQUENCE,BOOST_VMD_RETURN_ONLY_AFTER,BOOST_VMD_TYPE_LIST)));
  BOOST_TEST(!BOOST_VMD_IS_EMPTY(BOOST_VMD_ELEM(3,ASEQUENCE2,BOOST_VMD_RETURN_ONLY_AFTER,BOOST_VMD_TYPE_LIST)));
  BOOST_TEST(BOOST_VMD_IS_EMPTY(BOOST_VMD_ELEM(3,ASEQUENCE3,BOOST_VMD_RETURN_ONLY_AFTER,BOOST_VMD_TYPE_LIST)));
  BOOST_TEST(BOOST_VMD_IS_EMPTY(BOOST_VMD_ELEM(0,ASEQUENCE4,BOOST_VMD_RETURN_ONLY_AFTER,BOOST_VMD_TYPE_LIST)));
  BOOST_TEST(!BOOST_VMD_IS_EMPTY(BOOST_VMD_ELEM(0,ASEQUENCE5,BOOST_VMD_RETURN_ONLY_AFTER,BOOST_VMD_TYPE_LIST)));
  
  BOOST_TEST(BOOST_VMD_IS_EMPTY(BOOST_VMD_ELEM(0,anything,BOOST_VMD_RETURN_ONLY_AFTER,BOOST_VMD_TYPE_LIST)));
  BOOST_TEST_EQ(BOOST_VMD_ELEM(0,A_LIST_PLUS,BOOST_VMD_RETURN_ONLY_AFTER,BOOST_VMD_TYPE_LIST),456);
  BOOST_TEST(BOOST_VMD_IS_EMPTY(BOOST_VMD_ELEM(0,PLUS_ALIST,BOOST_VMD_RETURN_ONLY_AFTER,BOOST_VMD_TYPE_LIST)));
  BOOST_TEST(BOOST_VMD_IS_EMPTY(BOOST_VMD_ELEM(0,JDATA,BOOST_VMD_RETURN_ONLY_AFTER,BOOST_VMD_TYPE_LIST)));
  BOOST_TEST(!BOOST_VMD_IS_EMPTY(BOOST_VMD_ELEM(0,KDATA,BOOST_VMD_RETURN_ONLY_AFTER,BOOST_VMD_TYPE_LIST)));
  BOOST_TEST(!BOOST_VMD_IS_EMPTY(BOOST_VMD_ELEM(0,ELISTP,BOOST_VMD_RETURN_ONLY_AFTER,BOOST_VMD_TYPE_LIST)));
  BOOST_TEST_EQ(BOOST_VMD_ELEM(0,BOOST_PP_SEQ_ELEM(0,A_SEQ),BOOST_VMD_RETURN_ONLY_AFTER,BOOST_VMD_TYPE_LIST),3);
  
#else

BOOST_ERROR("No variadic macro support");
  
#endif

  return boost::report_errors();
  
  }
