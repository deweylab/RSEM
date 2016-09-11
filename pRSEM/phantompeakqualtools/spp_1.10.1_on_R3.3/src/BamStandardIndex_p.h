// ***************************************************************************
// BamStandardIndex.h (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 19 November 2010 (DB)
// ---------------------------------------------------------------------------
// Provides index operations for the standardized BAM index format (".bai")
// ***************************************************************************

#ifndef BAM_STANDARD_INDEX_FORMAT_H
#define BAM_STANDARD_INDEX_FORMAT_H

//  -------------
//  W A R N I N G
//  -------------
//
// This file is not part of the BamTools API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.

#include <BamAux.h>
#include <BamIndex.h>
#include <map>
#include <string>
#include <vector>

namespace BamTools {

class BamAlignment;

namespace Internal {

// BAM index constants
const int MAX_BIN        = 37450;    // =(8^6-1)/7+1
const int BAM_LIDX_SHIFT = 14;

// --------------------------------------------------
// BamStandardIndex data structures & typedefs
struct Chunk {

    // data members
    uint64_t Start;
    uint64_t Stop;

    // constructor
    Chunk(const uint64_t& start = 0,
	  const uint64_t& stop = 0)
	: Start(start)
	, Stop(stop)
    { }
};

inline
bool ChunkLessThan(const Chunk& lhs, const Chunk& rhs) {
    return lhs.Start < rhs.Start;
}

typedef std::vector<Chunk> ChunkVector;
typedef std::map<uint32_t, ChunkVector> BamBinMap;
typedef std::vector<uint64_t> LinearOffsetVector;

struct ReferenceIndex {

    // data members
    BamBinMap Bins;
    LinearOffsetVector Offsets;
    bool HasAlignments;

    // constructor
    ReferenceIndex(const BamBinMap& binMap           = BamBinMap(),
		   const LinearOffsetVector& offsets = LinearOffsetVector(),
		   const bool hasAlignments          = false)
	: Bins(binMap)
	, Offsets(offsets)
	, HasAlignments(hasAlignments)
    { }
};

typedef std::map<int32_t, ReferenceIndex> BamStandardIndexData;

class BamStandardIndex : public BamIndex {

    // ctor & dtor
    public:
	BamStandardIndex(BamTools::BgzfData* bgzf, BamTools::BamReader* reader);
	~BamStandardIndex(void);

    // interface (implements BamIndex virtual methods)
    public:
	// creates index data (in-memory) from current reader data
	bool Build(void);
	// returns supported file extension
	const std::string Extension(void) const { return std::string(".bai"); }
	// returns whether reference has alignments or no
	bool HasAlignments(const int& referenceID) const;
	// attempts to use index to jump to region; returns success/fail
	// a "successful" jump indicates no error, but not whether this region has data
	//   * thus, the method sets a flag to indicate whether there are alignments
	//     available after the jump position
	bool Jump(const BamTools::BamRegion& region, bool* hasAlignmentsInRegion);
    public:
	// clear all current index offset data in memory
	void ClearAllData(void);
	// return file position after header metadata
	const off_t DataBeginOffset(void) const;
	// return true if all index data is cached
	bool HasFullDataCache(void) const;
	// clears index data from all references except the first
	void KeepOnlyFirstReferenceOffsets(void);
	// load index data for all references, return true if loaded OK
	// @saveData - save data in memory if true, just read & discard if false
	bool LoadAllReferences(bool saveData = true);
	// load first reference from file, return true if loaded OK
	// @saveData - save data in memory if true, just read & discard if false
	bool LoadFirstReference(bool saveData = true);
	// load header data from index file, return true if loaded OK
	bool LoadHeader(void);
	// position file pointer to first reference begin, return true if skipped OK
	bool SkipToFirstReference(void);
	// write index reference data
	bool WriteAllReferences(void);
	// write index header data
	bool WriteHeader(void);

    // 'internal' methods
    public:

	// -----------------------
	// index file operations

	// check index file magic number, return true if OK
	bool CheckMagicNumber(void);
	// check index file version, return true if OK
	bool CheckVersion(void);
	// load a single index bin entry from file, return true if loaded OK
	// @saveData - save data in memory if true, just read & discard if false
	bool LoadBin(ReferenceIndex& refEntry, bool saveData = true);
	bool LoadBins(ReferenceIndex& refEntry, bool saveData = true);
	// load a single index bin entry from file, return true if loaded OK
	// @saveData - save data in memory if true, just read & discard if false
	bool LoadChunk(ChunkVector& chunks, bool saveData = true);
	bool LoadChunks(ChunkVector& chunks, bool saveData = true);
	// load a single index linear offset entry from file, return true if loaded OK
	// @saveData - save data in memory if true, just read & discard if false
	bool LoadLinearOffsets(ReferenceIndex& refEntry, bool saveData = true);
	// load a single reference from file, return true if loaded OK
	// @saveData - save data in memory if true, just read & discard if false
	bool LoadReference(const int& refId, bool saveData = true);
	// loads number of references, return true if loaded OK
	bool LoadReferenceCount(int& numReferences);
	// position file pointer to desired reference begin, return true if skipped OK
	bool SkipToReference(const int& refId);
	// write index data for bin to new index file
	bool WriteBin(const uint32_t& binId, const ChunkVector& chunks);
	// write index data for bins to new index file
	bool WriteBins(const BamBinMap& bins);
	// write index data for chunk entry to new index file
	bool WriteChunk(const Chunk& chunk);
	// write index data for chunk entry to new index file
	bool WriteChunks(const ChunkVector& chunks);
	// write index data for linear offsets entry to new index file
	bool WriteLinearOffsets(const LinearOffsetVector& offsets);
	// write index data single reference to new index file
	bool WriteReference(const ReferenceIndex& refEntry);

	// -----------------------
	// index data operations

	// calculate bins that overlap region
	int BinsFromRegion(const BamRegion& region,
			   const bool isRightBoundSpecified,
			   uint16_t bins[MAX_BIN]);
	// clear all index offset data for desired reference
	void ClearReferenceOffsets(const int& refId);
	// calculates offset(s) for a given region
	bool GetOffsets(const BamRegion& region,
			const bool isRightBoundSpecified,
			std::vector<int64_t>& offsets,
			bool* hasAlignmentsInRegion);
	// returns true if index cache has data for desired reference
	bool IsDataLoaded(const int& refId) const;
	// clears index data from all references except the one specified
	void KeepOnlyReferenceOffsets(const int& refId);
	// simplifies index by merging 'chunks'
	void MergeChunks(void);
	// saves BAM bin entry for index
	void SaveBinEntry(BamBinMap& binMap,
			  const uint32_t& saveBin,
			  const uint64_t& saveOffset,
			  const uint64_t& lastOffset);
	// saves linear offset entry for index
	void SaveLinearOffset(LinearOffsetVector& offsets,
			      const BamAlignment& bAlignment,
			      const uint64_t& lastOffset);
	// initializes index data structure to hold @count references
	void SetReferenceCount(const int& count);

    // data members
    private:

	BamStandardIndexData m_indexData;
	off_t m_dataBeginOffset;
	bool  m_hasFullDataCache;
	bool  m_isBigEndian;
};

} // namespace Internal
} // namespace BamTools

#endif // BAM_STANDARD_INDEX_FORMAT_H
