// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME XSTree_Dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "include/G2PRand.hh"
#include "include/HRSTransform_TCSNHCS.hh"
#include "XSModel/XSModel.hh"
#include "include/ACCInc.h"
#include "include/SHMSXSTree.h"
#include "include/HMSXSTree.h"
#include "include/ReadHMS.h"
#include "include/ReadSingleArm.h"
#include "include/ExtractAcceptance.h"
#include "include/ACCTools.h"
#include "include/ReadSHMS.h"
#include "include/XSTree.h"
#include "include/ACCIncBin.h"

// Header files passed via #pragma extra_include

namespace G2PRand {
   namespace ROOT {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *G2PRand_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("G2PRand", 0 /*version*/, "G2PRand.hh", 13,
                     ::ROOT::Internal::DefineBehavior((void*)0,(void*)0),
                     &G2PRand_Dictionary, 0);
         return &instance;
      }
      // Insure that the inline function is _not_ optimized away by the compiler
      ::ROOT::TGenericClassInfo *(*_R__UNIQUE_DICT_(InitFunctionKeeper))() = &GenerateInitInstance;  
      // Static variable to force the class initialization
      static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstance(); R__UseDummy(_R__UNIQUE_DICT_(Init));

      // Dictionary for non-ClassDef classes
      static TClass *G2PRand_Dictionary() {
         return GenerateInitInstance()->GetClass();
      }

   }
}

namespace Transform {
   namespace ROOT {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *Transform_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("Transform", 0 /*version*/, "HRSTransform_TCSNHCS.hh", 11,
                     ::ROOT::Internal::DefineBehavior((void*)0,(void*)0),
                     &Transform_Dictionary, 0);
         return &instance;
      }
      // Insure that the inline function is _not_ optimized away by the compiler
      ::ROOT::TGenericClassInfo *(*_R__UNIQUE_DICT_(InitFunctionKeeper))() = &GenerateInitInstance;  
      // Static variable to force the class initialization
      static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstance(); R__UseDummy(_R__UNIQUE_DICT_(Init));

      // Dictionary for non-ClassDef classes
      static TClass *Transform_Dictionary() {
         return GenerateInitInstance()->GetClass();
      }

   }
}

namespace ElasModel {
   namespace ROOT {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *ElasModel_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("ElasModel", 0 /*version*/, "XSModel.hh", 9,
                     ::ROOT::Internal::DefineBehavior((void*)0,(void*)0),
                     &ElasModel_Dictionary, 0);
         return &instance;
      }
      // Insure that the inline function is _not_ optimized away by the compiler
      ::ROOT::TGenericClassInfo *(*_R__UNIQUE_DICT_(InitFunctionKeeper))() = &GenerateInitInstance;  
      // Static variable to force the class initialization
      static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstance(); R__UseDummy(_R__UNIQUE_DICT_(Init));

      // Dictionary for non-ClassDef classes
      static TClass *ElasModel_Dictionary() {
         return GenerateInitInstance()->GetClass();
      }

   }
}

namespace PBosted {
   namespace ROOT {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *PBosted_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("PBosted", 0 /*version*/, "XSModel.hh", 38,
                     ::ROOT::Internal::DefineBehavior((void*)0,(void*)0),
                     &PBosted_Dictionary, 0);
         return &instance;
      }
      // Insure that the inline function is _not_ optimized away by the compiler
      ::ROOT::TGenericClassInfo *(*_R__UNIQUE_DICT_(InitFunctionKeeper))() = &GenerateInitInstance;  
      // Static variable to force the class initialization
      static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstance(); R__UseDummy(_R__UNIQUE_DICT_(Init));

      // Dictionary for non-ClassDef classes
      static TClass *PBosted_Dictionary() {
         return GenerateInitInstance()->GetClass();
      }

   }
}

namespace ROOT {
   static TClass *SHMSXSTree_Dictionary();
   static void SHMSXSTree_TClassManip(TClass*);
   static void delete_SHMSXSTree(void *p);
   static void deleteArray_SHMSXSTree(void *p);
   static void destruct_SHMSXSTree(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SHMSXSTree*)
   {
      ::SHMSXSTree *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::SHMSXSTree));
      static ::ROOT::TGenericClassInfo 
         instance("SHMSXSTree", "SHMSXSTree.h", 36,
                  typeid(::SHMSXSTree), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &SHMSXSTree_Dictionary, isa_proxy, 4,
                  sizeof(::SHMSXSTree) );
      instance.SetDelete(&delete_SHMSXSTree);
      instance.SetDeleteArray(&deleteArray_SHMSXSTree);
      instance.SetDestructor(&destruct_SHMSXSTree);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SHMSXSTree*)
   {
      return GenerateInitInstanceLocal((::SHMSXSTree*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SHMSXSTree*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *SHMSXSTree_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::SHMSXSTree*)0x0)->GetClass();
      SHMSXSTree_TClassManip(theClass);
   return theClass;
   }

   static void SHMSXSTree_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *HMSXSTree_Dictionary();
   static void HMSXSTree_TClassManip(TClass*);
   static void delete_HMSXSTree(void *p);
   static void deleteArray_HMSXSTree(void *p);
   static void destruct_HMSXSTree(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::HMSXSTree*)
   {
      ::HMSXSTree *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::HMSXSTree));
      static ::ROOT::TGenericClassInfo 
         instance("HMSXSTree", "HMSXSTree.h", 36,
                  typeid(::HMSXSTree), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &HMSXSTree_Dictionary, isa_proxy, 4,
                  sizeof(::HMSXSTree) );
      instance.SetDelete(&delete_HMSXSTree);
      instance.SetDeleteArray(&deleteArray_HMSXSTree);
      instance.SetDestructor(&destruct_HMSXSTree);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::HMSXSTree*)
   {
      return GenerateInitInstanceLocal((::HMSXSTree*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::HMSXSTree*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *HMSXSTree_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::HMSXSTree*)0x0)->GetClass();
      HMSXSTree_TClassManip(theClass);
   return theClass;
   }

   static void HMSXSTree_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *XSTree_Dictionary();
   static void XSTree_TClassManip(TClass*);
   static void delete_XSTree(void *p);
   static void deleteArray_XSTree(void *p);
   static void destruct_XSTree(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::XSTree*)
   {
      ::XSTree *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::XSTree));
      static ::ROOT::TGenericClassInfo 
         instance("XSTree", "XSTree.h", 36,
                  typeid(::XSTree), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &XSTree_Dictionary, isa_proxy, 4,
                  sizeof(::XSTree) );
      instance.SetDelete(&delete_XSTree);
      instance.SetDeleteArray(&deleteArray_XSTree);
      instance.SetDestructor(&destruct_XSTree);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::XSTree*)
   {
      return GenerateInitInstanceLocal((::XSTree*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::XSTree*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *XSTree_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::XSTree*)0x0)->GetClass();
      XSTree_TClassManip(theClass);
   return theClass;
   }

   static void XSTree_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrapper around operator delete
   static void delete_SHMSXSTree(void *p) {
      delete ((::SHMSXSTree*)p);
   }
   static void deleteArray_SHMSXSTree(void *p) {
      delete [] ((::SHMSXSTree*)p);
   }
   static void destruct_SHMSXSTree(void *p) {
      typedef ::SHMSXSTree current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SHMSXSTree

namespace ROOT {
   // Wrapper around operator delete
   static void delete_HMSXSTree(void *p) {
      delete ((::HMSXSTree*)p);
   }
   static void deleteArray_HMSXSTree(void *p) {
      delete [] ((::HMSXSTree*)p);
   }
   static void destruct_HMSXSTree(void *p) {
      typedef ::HMSXSTree current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::HMSXSTree

namespace ROOT {
   // Wrapper around operator delete
   static void delete_XSTree(void *p) {
      delete ((::XSTree*)p);
   }
   static void deleteArray_XSTree(void *p) {
      delete [] ((::XSTree*)p);
   }
   static void destruct_XSTree(void *p) {
      typedef ::XSTree current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::XSTree

namespace {
  void TriggerDictionaryInitialization_XSTree_Dict_Impl() {
    static const char* headers[] = {
"include/G2PRand.hh",
"include/HRSTransform_TCSNHCS.hh",
"XSModel/XSModel.hh",
"include/ACCInc.h",
"include/SHMSXSTree.h",
"include/HMSXSTree.h",
"include/ReadHMS.h",
"include/ReadSingleArm.h",
"include/ExtractAcceptance.h",
"include/ACCTools.h",
"include/ReadSHMS.h",
"include/XSTree.h",
"include/ACCIncBin.h",
0
    };
    static const char* includePaths[] = {
"XSModel",
"XSModel/include",
"/site/12gev_phys/2.1/Linux_CentOS7.7.1908-x86_64-gcc4.8.5/root/6.10.02/include",
"/u/site/12gev_phys/2.1/Linux_CentOS7.2.1511-x86_64-gcc4.8.5/root/6.10.02/include",
"/u/group/c-polhe3/Users/murchhana/CreateXSTree/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "XSTree_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$include/SHMSXSTree.h")))  SHMSXSTree;
class __attribute__((annotate("$clingAutoload$include/HMSXSTree.h")))  HMSXSTree;
class __attribute__((annotate("$clingAutoload$include/XSTree.h")))  XSTree;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "XSTree_Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "include/G2PRand.hh"
#include "include/HRSTransform_TCSNHCS.hh"
#include "XSModel/XSModel.hh"
#include "include/ACCInc.h"
#include "include/SHMSXSTree.h"
#include "include/HMSXSTree.h"
#include "include/ReadHMS.h"
#include "include/ReadSingleArm.h"
#include "include/ExtractAcceptance.h"
#include "include/ACCTools.h"
#include "include/ReadSHMS.h"
#include "include/XSTree.h"
#include "include/ACCIncBin.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"HMSXSTree", payloadCode, "@",
"SHMSXSTree", payloadCode, "@",
"XSTree", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("XSTree_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_XSTree_Dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_XSTree_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_XSTree_Dict() {
  TriggerDictionaryInitialization_XSTree_Dict_Impl();
}
