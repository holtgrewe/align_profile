// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
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

#ifndef SANDBOX_ALIGN_PROFILE_INCLUDE_SEQAN_ALIGN_PROFILE_ALIGN_PROFILE_SCORE_PROFILE_SEQ_H_
#define SANDBOX_ALIGN_PROFILE_INCLUDE_SEQAN_ALIGN_PROFILE_ALIGN_PROFILE_SCORE_PROFILE_SEQ_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

struct ProfileSeqScore_;
typedef Tag<ProfileSeqScore_> ProfileSeqScore;

template <typename TValue>
class Score<TValue, ProfileSeqScore>;

template <typename TValue, typename TString>
inline void
assignProfile(Score<TValue, ProfileSeqScore> & me, TString const & profile);

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class ProfileSeq Score
// ----------------------------------------------------------------------------

struct ProfileSeqScore_;
typedef Tag<ProfileSeqScore_> ProfileSeqScore;

/**
.Spec.ProfileSeq Score
..summary:Score for profile-to-sequence alignments.
*/

template <typename TValue>
class Score<TValue, ProfileSeqScore>
{
public:
	String<TValue> consensusSet;		// Is the alphabet character part of the consensus set for each column

    Score() {}

    // Construct given a profile string.
    template <typename TProfile>
    explicit
    Score(TProfile const & profile)
    {
        assignProfile(*this, profile);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// --------------------------------------------------------------------------
// Metafunction SequenceEntryForScore                      [ProfileSeq Score]
// --------------------------------------------------------------------------

// Returns the type that holds a sequence entry.  This is used for abstracting away the access to sequence characters.

template <typename TValue, typename TSequence>
struct SequenceEntryForScore<Score<TValue, ProfileSeqScore>, TSequence>
{
    typedef ConsensusScoreSequenceEntry<TSequence> Type;
};

template <typename TValue, typename TSequence>
struct SequenceEntryForScore<Score<TValue, ProfileSeqScore> const, TSequence> :
       SequenceEntryForScore<Score<TValue, ProfileSeqScore>, TSequence>
{};

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function sequenceEntryForScore()                        [ProfileSeq Score]
// --------------------------------------------------------------------------

template <typename TScoreValue, typename TSequence, typename TPosition>
inline ConsensusScoreSequenceEntry<TSequence>
sequenceEntryForScore(Score<TScoreValue, ProfileSeqScore> const & /*sScheme*/,
                      TSequence const & seq, TPosition pos)
{
    return ConsensusScoreSequenceEntry<TSequence>(seq, pos);
}

// --------------------------------------------------------------------------
// Function assignProfile()                                [ProfileSeq Score]
// --------------------------------------------------------------------------

template <typename TValue, typename TString>
inline void
assignProfile(Score<TValue, ProfileSeqScore> & me,
			  TString const & profile)
{
	typedef typename Size<TString>::Type TSize;
	TSize alphSize = ValueSize<typename Value<TString>::Type>::VALUE;
	resize(me.consensusSet, alphSize * length(profile));

	typedef typename Iterator<TString const, Standard>::Type TIter;
	typedef typename Iterator<String<TValue>, Standard>::Type TConsSetIter;
	TConsSetIter itConsSet = begin(me.consensusSet, Standard());
	TIter it = begin(profile, Standard());
	TIter itEnd = end(profile, Standard());
	TSize maxCount = 0;
	for(;it!=itEnd;++it) {
		maxCount = 0;
		for(TSize i = 0; i<alphSize; ++i)
			if ((TSize)(*it).count[i] > maxCount)
                maxCount = (*it).count[i];
		for(TSize i = 0; i<alphSize; ++i, ++itConsSet)
			*itConsSet = ((TSize)(*it).count[i] == maxCount)? 0 : (-SEQAN_CONSENSUS_UNITY);
	}
}

// --------------------------------------------------------------------------
// Function scoreGapExtendHorizontal()                     [ProfileSeq Score]
// --------------------------------------------------------------------------

template <typename TValue, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendHorizontal(
	Score<TValue, ProfileSeqScore> const & me,
    ConsensusScoreSequenceEntry<TSeq1> const & entry1,
    ConsensusScoreSequenceEntry<TSeq2> const & entry2)
{
	return ((int)position(entry2) < 0) ? -SEQAN_CONSENSUS_UNITY : me.consensusSet[position(entry1) * (ValueSize<typename Value<TSeq1>::Type>::VALUE) + (ValueSize<typename Value<TSeq1>::Type>::VALUE - 1)];
}

// --------------------------------------------------------------------------
// Function scoreGapOpenHorizontal()                       [ProfileSeq Score]
// --------------------------------------------------------------------------

template <typename TValue, typename TSeq1, typename TSeq2>
inline TValue
scoreGapOpenHorizontal(
	Score<TValue, ProfileSeqScore> const & me,
    ConsensusScoreSequenceEntry<TSeq1> const & entry1,
    ConsensusScoreSequenceEntry<TSeq2> const & entry2)
{
	return ((int)position(entry2) < 0) ? -2 * SEQAN_CONSENSUS_UNITY : 2 * me.consensusSet[position(entry1) * (ValueSize<typename Value<TSeq1>::Type>::VALUE) + (ValueSize<typename Value<TSeq1>::Type>::VALUE - 1)];
}

// --------------------------------------------------------------------------
// Function scoreGapExtendVertical()                       [ProfileSeq Score]
// --------------------------------------------------------------------------

template <typename TValue, typename TSeq1, typename TSeq2>
inline TValue
scoreGapOpenVertical(
	Score<TValue, ProfileSeqScore> const &,
    ConsensusScoreSequenceEntry<TSeq1> const & /*entry1*/,
    ConsensusScoreSequenceEntry<TSeq2> const & /*entry2*/)
{
	return -2 * SEQAN_CONSENSUS_UNITY;
}

// --------------------------------------------------------------------------
// Function scoreGapOpenVertical()                         [ProfileSeq Score]
// --------------------------------------------------------------------------

template <typename TValue, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendVertical(
	Score<TValue, ProfileSeqScore> const &,
    ConsensusScoreSequenceEntry<TSeq1> const & /*entry1*/,
    ConsensusScoreSequenceEntry<TSeq2> const & /*entry2*/)
{
	return -SEQAN_CONSENSUS_UNITY;
}

// --------------------------------------------------------------------------
// Function score()                                        [ProfileSeq Score]
// --------------------------------------------------------------------------

template <typename TValue, typename TSeq1, typename TSeq2>
inline TValue
score(Score<TValue, ProfileSeqScore> const & me,
      ConsensusScoreSequenceEntry<TSeq1> const & entry1,
      ConsensusScoreSequenceEntry<TSeq2> const & entry2)
{
	return me.consensusSet[position(entry1) * (ValueSize<typename Value<TSeq1>::Type>::VALUE) + ordValue(value(entry2))];
}

}  // namespace seqan

#endif  // #ifndef SANDBOX_ALIGN_PROFILE_INCLUDE_SEQAN_ALIGN_PROFILE_ALIGN_PROFILE_SCORE_PROFILE_SEQ_H_
