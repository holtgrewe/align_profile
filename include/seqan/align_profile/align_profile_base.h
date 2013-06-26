// ==========================================================================
//                               align_profile
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef SANDBOX_ALIGN_PROFILE_INCLUDE_SEQAN_ALIGN_PROFILE_ALIGN_PROFILE_BASE_H_
#define SANDBOX_ALIGN_PROFILE_INCLUDE_SEQAN_ALIGN_PROFILE_ALIGN_PROFILE_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Concept:MyConcept
..
..include:seqan/align_profile.h
 */

// If your define a lot of very generic functions for our concept, consider
// putting it into its own directory.

/**
.Class.MyClass
..concept:Concept.MyConcept
..summary:This is my class.
..cat:My Classes
..signature:MyClass<TSpec>
..param.TSpec:Tag to select the specialization.
...default:MyTag
..include:seqan/align_profile.h
 */

// NOTE: Assigning classes to concepts is optional.

// struct Our_;
// typedef Tag<Our_> Our;

// template <typename TSpec = Our>
// class MyClass;

/*
.Spec.Our MyClass
..cat:My Classes
..general:Class.MyClass
..summary:This is the "our" specialization of my class!
..signature:MyClass<Our>
..include:seqan/align_profile.h
 
.Memfunc.Our MyClass#MyClass
..cat:My Classes
..class:Spec.Our MyClass
..signature:Class()
..signature:Class(foo, barBaz)
..param.foo:A foo parameter.
...type:Spec.CharString
..param.barBaz:Another parameter.
...type:nolink:$int$
 */
	
// template <>
// class MyClass<Our>
// {
// public:
//     // ...
// };

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

///.Metafunction.Value.param.T.value:Class.MyClass

// template <typename TSpec>
// struct Value<MyClass<TSpec> >
// {
//     typedef int Type;
// };

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function someFunction()
// ----------------------------------------------------------------------------

/**
.Function.someFunction
..concept:Concept.MyConcept
..cat:My Classes
..signature:someFunction(obj)
..summary:Executes some functionality on a @Class.MyClass@.
..param.obj:Object to call function on.
...type:Class.MyClass
..returns:$void$
..include:seqan/align_profile.h
 */

// NOTE: Functions can belong to concepts but do not have to.  You can use
//       function documentations and link them to concepts in a concept_name.h
//       header to document your concept.

// template <typename TSpec>
// inline void
// someFunction(MyClass<TSpec> const & /*obj*/)
// {
//     return;
// }

}  // namespace seqan

#endif  // SANDBOX_ALIGN_PROFILE_INCLUDE_SEQAN_ALIGN_PROFILE_ALIGN_PROFILE_BASE_H_
