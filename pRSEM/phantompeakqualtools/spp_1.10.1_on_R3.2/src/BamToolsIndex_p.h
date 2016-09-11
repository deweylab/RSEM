// ***************************************************************************
// BamToolsIndex.h (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 19 November 2010 (DB)
// ---------------------------------------------------------------------------
// Provides index operations for the BamTools index format (".bti")
// ***************************************************************************

#ifndef BAMTOOLS_INDEX_FORMAT_H
#define BAMTOOLS_INDEX_FORMAT_H

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

namespace Internal {

// individual index offset entry
struct BamToolsIndexEntry {

    // data members
    int32_t MaxEndPosition;
    int64_t StartOffset;
    int32_t StartPosition;

    // ctor
    BamToolsIndexEntry(const int32_t& maxEndPosition = 0,
		       const int64_t& startOffset    = 0,
		       const int32_t& startPosition  = 0)
	: MaxEndPosition(maxEndPosition)
	, StartOffset(startOffset)
	, StartPosition(startPosition)
    { }
};

// reference index entry
struct BamToolsReferenceEntry {

    // data members
    bool HasAlignments;
    std::vector<BamToolsIndexEntry> Offsets;

    // ctor
    BamToolsReferenceEntry(void)
	: HasAlignments(false)
    { }
};

// the actual index data structure
typedef std::map<int, BamToolsReferenceEntry> BamToolsIndexData;

class BamToolsIndex : public BamIndex {

    // keep a list of any supported versions here
    // (might be useful later to handle any 'legacy' versions if the format changes)
    // listed for example like: BTI_1_0 = 1, BTI_1_1 = 2, BTI_1_2 = 3, BTI_2_0 = 4, and so on
    //
    // so a change introduced in (hypothetical) BTI_1_2 would be handled from then on by:
    //
    // if ( indexVersion >= BTI_1_2 )
    //   do something new
    // else
    //   do the old thing
    enum Version { BTI_1_0 = 1
		 , BTI_1_1
		 , BTI_1_2
		 };


    // ctor & dtor
    public:
	BamToolsIndex(BamTools::BgzfData* bgzf, BamTools::BamReader* reader);
	~BamToolsIndex(void);

    // interface (implements BamIndex virtual methods)
    public:
	// creates index data (in-memory) from current reader data
	bool Build(void);
	// returns supported file extension
	const std::string Extension(void) const { return std::string(".bti"); }
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
	// return true if FILE* is open
	bool IsOpen(void) const;
	// load a single index entry from file, return true if loaded OK
	// @saveData - save data in memory if true, just read & discard if false
	bool LoadIndexEntry(const int& refId, bool saveData = true);
	// load a single reference from file, return true if loaded OK
	// @saveData - save data in memory if true, just read & discard if false
	bool LoadReference(const int& refId, bool saveData = true);
	// loads number of references, return true if loaded OK
	bool LoadReferenceCount(int& numReferences);
	// position file pointer to desired reference begin, return true if skipped OK
	bool SkipToReference(const int& refId);
	// write current reference index data to new index file
	bool WriteReferenceEntry(const BamToolsReferenceEntry& refEntry);
	// write current index offset entry to new index file
	bool WriteIndexEntry(const BamToolsIndexEntry& entry);

	// -----------------------
	// index data operations

	// clear all index offset data for desired reference
	void ClearReferenceOffsets(const int& refId);
	// calculate BAM file offset for desired region
	// return true if no error (*NOT* equivalent to "has alignments or valid offset")
	//   check @hasAlignmentsInRegion to determine this status
	// @region - target region
	// @offset - resulting seek target
	// @hasAlignmentsInRegion - sometimes a file just lacks data in region, this flag indicates that status
	bool GetOffset(const BamRegion& region, int64_t& offset, bool* hasAlignmentsInRegion);
	// returns true if index cache has data for desired reference
	bool IsDataLoaded(const int& refId) const;
	// clears index data from all references except the one specified
	void KeepOnlyReferenceOffsets(const int& refId);
	// saves an index offset entry in memory
	void SaveOffsetEntry(const int& refId, const BamToolsIndexEntry& entry);
	// pre-allocates size for offset vector
	void SetOffsetCount(const int& refId, const int& offsetCount);
	// initializes index data structure to hold @count references
	void SetReferenceCount(const int& count);

    // data members
    private:
	int32_t           m_blockSize;
	BamToolsIndexData m_indexData;
	off_t             m_dataBeginOffset;
	bool              m_hasFullDataCache;
	bool              m_isBigEndian;
	int32_t           m_inputVersion; // Version is serialized as int
	Version           m_outputVersion;
};

} // namespace Internal
} // namespace BamTools

#endif // BAMTOOLS_INDEX_FORMAT_H
