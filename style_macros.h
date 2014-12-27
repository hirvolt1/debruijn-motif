// Copyright 2012 Kalle Karhu FREE AS A FLY
#ifndef GTEST_STYLE_MACROS_H_
#define GTEST_STYLE_MACROS_H_
// A macro to disallow the copy constructor and operator= functions
// This should be used in the private: declarations for class
// Source: Google C++ Style Guide
#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
  TypeName(const TypeName&);               \
  void operator=(const TypeName&)

#endif  // GTEST_STYLE_MACROS_H_
