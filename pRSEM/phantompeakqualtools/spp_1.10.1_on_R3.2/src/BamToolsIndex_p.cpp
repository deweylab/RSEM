// ***************************************************************************
// BamToolsIndex.cpp (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 22 November 2010 (DB)
// ---------------------------------------------------------------------------
// Provides index operations for the BamTools index format (".bti")
// ***************************************************************************

#include <BamAlignment.h>
#include <BamReader.h>
#include <BGZF.h>
#include <BamToolsIndex_p.h>
using namespace BamTools;
using namespace BamTools::Internal;

#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <map>
using namespace std;

BamToolsIndex::BamToolsIndex(BgzfData* bgzf, BamReader* reader)
    : BamIndex(bgzf, reader)
    , m_blockSize(1000)
    , m_dataBeginOffset(0)
    , m_hasFullDataCache(false)
    , m_inputVersion(0)
    , m_outputVersion(BTI_1_2) // latest version - used for writing new index files
{
    m_isBigEndian = BamTools::SystemIsBigEndian();
}

// dtor
BamToolsIndex::~BamToolsIndex(void) {
    ClearAllData();
}

// creates index data (in-memory) from current reader data
bool BamToolsIndex::Build(void) {

    // be sure reader & BGZF file are valid & open for reading
    if ( m_reader == 0 || m_BGZF == 0 || !m_BGZF->IsOpen )
	return false;

    // move file pointer to beginning of alignments
    if ( !m_reader->Rewind() ) return false;

    // initialize index data structure with space for all references
    const int numReferences = (int)m_references.size();
    m_indexData.clear();
    m_hasFullDataCache = false;
    SetReferenceCount(numReferences);

    // set up counters and markers
    int32_t currentBlockCount      = 0;
    int64_t currentAlignmentOffset = m_BGZF->Tell();
    int32_t blockRefId             = 0;
    int32_t blockMaxEndPosition    = 0;
    int64_t blockStartOffset       = currentAlignmentOffset;
    int32_t blockStartPosition     = -1;

    // plow through alignments, storing index entries
    BamAlignment al;
    while ( m_reader->GetNextAlignmentCore(al) ) {

	// if block contains data (not the first time through) AND alignment is on a new reference
	if ( currentBlockCount > 0 && al.RefID != blockRefId ) {

	    // store previous data
	    BamToolsIndexEntry entry(blockMaxEndPosition, blockStartOffset, blockStartPosition);
	    SaveOffsetEntry(blockRefId, entry);

	    // intialize new block for current alignment's reference
	    currentBlockCount   = 0;
	    blockMaxEndPosition = al.GetEndPosition();
	    blockStartOffset    = currentAlignmentOffset;
	}

	// if beginning of block, save first alignment's refID & position
	if ( currentBlockCount == 0 ) {
	    blockRefId         = al.RefID;
	    blockStartPosition = al.Position;
	}

	// increment block counter
	++currentBlockCount;

	// check end position
	int32_t alignmentEndPosition = al.GetEndPosition();
	if ( alignmentEndPosition > blockMaxEndPosition )
	    blockMaxEndPosition = alignmentEndPosition;

	// if block is full, get offset for next block, reset currentBlockCount
	if ( currentBlockCount == m_blockSize ) {
	    BamToolsIndexEntry entry(blockMaxEndPosition, blockStartOffset, blockStartPosition);
	    SaveOffsetEntry(blockRefId, entry);
	    blockStartOffset  = m_BGZF->Tell();
	    currentBlockCount = 0;
	}

	// not the best name, but for the next iteration, this value will be the offset of the *current* alignment
	// necessary because we won't know if this next alignment is on a new reference until we actually read it
	currentAlignmentOffset = m_BGZF->Tell();
    }

    // store final block with data
    BamToolsIndexEntry entry(blockMaxEndPosition, blockStartOffset, blockStartPosition);
    SaveOffsetEntry(blockRefId, entry);

    // set flag
    m_hasFullDataCache = true;

    // return success/failure of rewind
    return m_reader->Rewind();
}

// check index file magic number, return true if OK
bool BamToolsIndex::CheckMagicNumber(void) {

    // see if index is valid BAM index
    char magic[4];
    size_t elementsRead = fread(magic, 1, 4, m_indexStream);
    if ( elementsRead != 4 ) return false;
    if ( strncmp(magic, "BTI\1", 4) != 0 ) {
	fprintf(stderr, "Problem with index file - invalid format.\n");
	return false;
    }

    // otherwise ok
    return true;
}

// check index file version, return true if OK
bool BamToolsIndex::CheckVersion(void) {

    // read version from file
    size_t elementsRead = fread(&m_inputVersion, sizeof(m_inputVersion), 1, m_indexStream);
    if ( elementsRead != 1 ) return false;
    if ( m_isBigEndian ) SwapEndian_32(m_inputVersion);

    // if version is negative, or zero
    if ( m_inputVersion <= 0 ) {
	fprintf(stderr, "Problem with index file - invalid version.\n");
	return false;
    }

    // if version is newer than can be supported by this version of bamtools
    else if ( m_inputVersion > m_outputVersion ) {
	fprintf(stderr, "Problem with index file - attempting to use an outdated version of BamTools with a newer index file.\n");
	fprintf(stderr, "Please update BamTools to a more recent version to support this index file.\n");
	return false;
    }

    // ------------------------------------------------------------------
    // check for deprecated, unsupported versions
    // (typically whose format did not accomodate a particular bug fix)

    else if ( (Version)m_inputVersion == BTI_1_0 ) {
	fprintf(stderr, "\nProblem with index file - this version of the index contains a bug related to accessing data near reference ends.\n");
	fprintf(stderr, "\nPlease run \'bamtools index -bti -in yourData.bam\' to generate an up-to-date BamToolsIndex.\n\n");
	return false;
    }

    else if ( (Version)m_inputVersion == BTI_1_1 ) {
	fprintf(stderr, "\nProblem with index file - this version of the index contains a bug related to handling empty references.\n");
	fprintf(stderr, "\nPlease run \'bamtools index -bti -in yourData.bam\' to generate an up-to-date BamToolsIndex.\n\n");
	return false;
    }

    // otherwise ok
    else return true;
}

// clear all current index offset data in memory
void BamToolsIndex::ClearAllData(void) {
    BamToolsIndexData::const_iterator indexIter = m_indexData.begin();
    BamToolsIndexData::const_iterator indexEnd  = m_indexData.end();
    for ( ; indexIter != indexEnd; ++indexIter ) {
	const int& refId = (*indexIter).first;
	ClearReferenceOffsets(refId);
    }
}

// clear all index offset data for desired reference
void BamToolsIndex::ClearReferenceOffsets(const int& refId) {
    if ( m_indexData.find(refId) == m_indexData.end() ) return;
    vector<BamToolsIndexEntry>& offsets = m_indexData[refId].Offsets;
    offsets.clear();
    m_hasFullDataCache = false;
}

// return file position after header metadata
const off_t BamToolsIndex::DataBeginOffset(void) const {
    return m_dataBeginOffset;
}

// calculate BAM file offset for desired region
// return true if no error (*NOT* equivalent to "has alignments or valid offset")
//   check @hasAlignmentsInRegion to determine this status
// @region - target region
// @offset - resulting seek target
// @hasAlignmentsInRegion - sometimes a file just lacks data in region, this flag indicates that status
// N.B. - ignores isRightBoundSpecified
bool BamToolsIndex::GetOffset(const BamRegion& region, int64_t& offset, bool* hasAlignmentsInRegion) {

    // return false if leftBound refID is not found in index data
    BamToolsIndexData::const_iterator indexIter = m_indexData.find(region.LeftRefID);
    if ( indexIter == m_indexData.end()) return false;

    // load index data for region if not already cached
    if ( !IsDataLoaded(region.LeftRefID) ) {
	bool loadedOk = true;
	loadedOk &= SkipToReference(region.LeftRefID);
	loadedOk &= LoadReference(region.LeftRefID);
	if ( !loadedOk ) return false;
    }

    // localize index data for this reference (& sanity check that data actually exists)
    indexIter = m_indexData.find(region.LeftRefID);
    if ( indexIter == m_indexData.end()) return false;
    const vector<BamToolsIndexEntry>& referenceOffsets = (*indexIter).second.Offsets;
    if ( referenceOffsets.empty() ) return false;

    // -------------------------------------------------------
    // calculate nearest index to jump to

    // save first offset
    offset = (*referenceOffsets.begin()).StartOffset;

    // iterate over offsets entries on this reference
    vector<BamToolsIndexEntry>::const_iterator offsetIter = referenceOffsets.begin();
    vector<BamToolsIndexEntry>::const_iterator offsetEnd  = referenceOffsets.end();
    for ( ; offsetIter != offsetEnd; ++offsetIter ) {
	const BamToolsIndexEntry& entry = (*offsetIter);
	// break if alignment 'entry' overlaps region
	if ( entry.MaxEndPosition >= region.LeftPosition ) break;
	offset = (*offsetIter).StartOffset;
    }

    // set flag based on whether an index entry was found for this region
    *hasAlignmentsInRegion = ( offsetIter != offsetEnd );

    // if cache mode set to none, dump the data we just loaded
    if (m_cacheMode == BamIndex::NoIndexCaching )
	ClearReferenceOffsets(region.LeftRefID);

    // return success
    return true;
}

// returns whether reference has alignments or no
bool BamToolsIndex::HasAlignments(const int& refId) const {

    BamToolsIndexData::const_iterator indexIter = m_indexData.find(refId);
    if ( indexIter == m_indexData.end()) return false;
    const BamToolsReferenceEntry& refEntry = (*indexIter).second;
    return refEntry.HasAlignments;
}

// return true if all index data is cached
bool BamToolsIndex::HasFullDataCache(void) const {
    return m_hasFullDataCache;
}

// returns true if index cache has data for desired reference
bool BamToolsIndex::IsDataLoaded(const int& refId) const {

    BamToolsIndexData::const_iterator indexIter = m_indexData.find(refId);
    if ( indexIter == m_indexData.end()) return false;
    const BamToolsReferenceEntry& refEntry = (*indexIter).second;

    if ( !refEntry.HasAlignments ) return true; // no data period

    // return whether offsets list contains data
    return !refEntry.Offsets.empty();
}

// attempts to use index to jump to region; returns success/fail
bool BamToolsIndex::Jump(const BamRegion& region, bool* hasAlignmentsInRegion) {

    // clear flag
    *hasAlignmentsInRegion = false;

    // check valid BamReader state
    if ( m_reader == 0 || m_BGZF == 0 || !m_reader->IsOpen() ) {
	fprintf(stderr, "ERROR: Could not jump: invalid BamReader state.\n");
	return false;
    }

    // make sure left-bound position is valid
    if ( region.LeftPosition > m_references.at(region.LeftRefID).RefLength )
	return false;

    // calculate nearest offset to jump to
    int64_t offset;
    if ( !GetOffset(region, offset, hasAlignmentsInRegion) ) {
	fprintf(stderr, "ERROR: Could not jump - unable to calculate offset for specified region.\n");
	return false;
    }

    // return success/failure of seek
    return m_BGZF->Seek(offset);
}

// clears index data from all references except the first
void BamToolsIndex::KeepOnlyFirstReferenceOffsets(void) {
    BamToolsIndexData::const_iterator indexBegin = m_indexData.begin();
    KeepOnlyReferenceOffsets( (*indexBegin).first );
}

// clears index data from all references except the one specified
void BamToolsIndex::KeepOnlyReferenceOffsets(const int& refId) {
    BamToolsIndexData::iterator mapIter = m_indexData.begin();
    BamToolsIndexData::iterator mapEnd  = m_indexData.end();
    for ( ; mapIter != mapEnd; ++mapIter ) {
	const int entryRefId = (*mapIter).first;
	if ( entryRefId != refId )
	    ClearReferenceOffsets(entryRefId);
    }
}

// load index data for all references, return true if loaded OK
bool BamToolsIndex::LoadAllReferences(bool saveData) {

    // skip if data already loaded
    if ( m_hasFullDataCache ) return true;

    // read in number of references
    int32_t numReferences;
    if ( !LoadReferenceCount(numReferences) ) return false;
    //SetReferenceCount(numReferences);

    // iterate over reference entries
    bool loadedOk = true;
    for ( int i = 0; i < numReferences; ++i )
	loadedOk &= LoadReference(i, saveData);

    // set flag
    if ( loadedOk && saveData )
	m_hasFullDataCache = true;

    // return success/failure of load
    return loadedOk;
}

// load header data from index file, return true if loaded OK
bool BamToolsIndex::LoadHeader(void) {

    // check magic number
    if ( !CheckMagicNumber() ) return false;

    // check BTI version
    if ( !CheckVersion() ) return false;

    // read in block size
    size_t elementsRead = fread(&m_blockSize, sizeof(m_blockSize), 1, m_indexStream);
    if ( elementsRead != 1 ) return false;
    if ( m_isBigEndian ) SwapEndian_32(m_blockSize);

    // store offset of beginning of data
    m_dataBeginOffset = ftell64(m_indexStream);

    // return success/failure of load
    return (elementsRead == 1);
}

// load a single index entry from file, return true if loaded OK
// @saveData - save data in memory if true, just read & discard if false
bool BamToolsIndex::LoadIndexEntry(const int& refId, bool saveData) {

    // read in index entry data members
    size_t elementsRead = 0;
    BamToolsIndexEntry entry;
    elementsRead += fread(&entry.MaxEndPosition, sizeof(entry.MaxEndPosition), 1, m_indexStream);
    elementsRead += fread(&entry.StartOffset,    sizeof(entry.StartOffset),    1, m_indexStream);
    elementsRead += fread(&entry.StartPosition,  sizeof(entry.StartPosition),  1, m_indexStream);
    if ( elementsRead != 3 ) {
	cerr << "Error reading index entry. Expected 3 elements, read in: " << elementsRead << endl;
	return false;
    }

    // swap endian-ness if necessary
    if ( m_isBigEndian ) {
	SwapEndian_32(entry.MaxEndPosition);
	SwapEndian_64(entry.StartOffset);
	SwapEndian_32(entry.StartPosition);
    }

    // save data
    if ( saveData )
	SaveOffsetEntry(refId, entry);

    // return success/failure of load
    return true;
}

// load a single reference from file, return true if loaded OK
// @saveData - save data in memory if true, just read & discard if false
bool BamToolsIndex::LoadFirstReference(bool saveData) {
    BamToolsIndexData::const_iterator indexBegin = m_indexData.begin();
    return LoadReference( (*indexBegin).first, saveData );
}

// load a single reference from file, return true if loaded OK
// @saveData - save data in memory if true, just read & discard if false
bool BamToolsIndex::LoadReference(const int& refId, bool saveData) {

    // read in number of offsets for this reference
    uint32_t numOffsets;
    size_t elementsRead = fread(&numOffsets, sizeof(numOffsets), 1, m_indexStream);
    if ( elementsRead != 1 ) return false;
    if ( m_isBigEndian ) SwapEndian_32(numOffsets);

    // initialize offsets container for this reference
    SetOffsetCount(refId, (int)numOffsets);

    // iterate over offset entries
    for ( unsigned int j = 0; j < numOffsets; ++j )
	LoadIndexEntry(refId, saveData);

    // return success/failure of load
    return true;
}

// loads number of references, return true if loaded OK
bool BamToolsIndex::LoadReferenceCount(int& numReferences) {

    size_t elementsRead = 0;

    // read reference count
    elementsRead += fread(&numReferences, sizeof(numReferences), 1, m_indexStream);
    if ( m_isBigEndian ) SwapEndian_32(numReferences);

    // return success/failure of load
    return ( elementsRead == 1 );
}

// saves an index offset entry in memory
void BamToolsIndex::SaveOffsetEntry(const int& refId, const BamToolsIndexEntry& entry) {
    BamToolsReferenceEntry& refEntry = m_indexData[refId];
    refEntry.HasAlignments = true;
    refEntry.Offsets.push_back(entry);
}

// pre-allocates size for offset vector
void BamToolsIndex::SetOffsetCount(const int& refId, const int& offsetCount) {
    BamToolsReferenceEntry& refEntry = m_indexData[refId];
    refEntry.Offsets.reserve(offsetCount);
    refEntry.HasAlignments = ( offsetCount > 0);
}

// initializes index data structure to hold @count references
void BamToolsIndex::SetReferenceCount(const int& count) {
    for ( int i = 0; i < count; ++i )
	m_indexData[i].HasAlignments = false;
}

// position file pointer to first reference begin, return true if skipped OK
bool BamToolsIndex::SkipToFirstReference(void) {
    BamToolsIndexData::const_iterator indexBegin = m_indexData.begin();
    return SkipToReference( (*indexBegin).first );
}

// position file pointer to desired reference begin, return true if skipped OK
bool BamToolsIndex::SkipToReference(const int& refId) {

    // attempt rewind
    if ( !Rewind() ) return false;

    // read in number of references
    int32_t numReferences;
    size_t elementsRead = fread(&numReferences, sizeof(numReferences), 1, m_indexStream);
    if ( elementsRead != 1 ) return false;
    if ( m_isBigEndian ) SwapEndian_32(numReferences);

    // iterate over reference entries
    bool skippedOk = true;
    int currentRefId = 0;
    while (currentRefId != refId) {
	skippedOk &= LoadReference(currentRefId, false);
	++currentRefId;
    }

    // return success/failure of skip
    return skippedOk;
}

// write header to new index file
bool BamToolsIndex::WriteHeader(void) {

    size_t elementsWritten = 0;

    // write BTI index format 'magic number'
    elementsWritten += fwrite("BTI\1", 1, 4, m_indexStream);

    // write BTI index format version
    int32_t currentVersion = (int32_t)m_outputVersion;
    if ( m_isBigEndian ) SwapEndian_32(currentVersion);
    elementsWritten += fwrite(&currentVersion, sizeof(currentVersion), 1, m_indexStream);

    // write block size
    int32_t blockSize = m_blockSize;
    if ( m_isBigEndian ) SwapEndian_32(blockSize);
    elementsWritten += fwrite(&blockSize, sizeof(blockSize), 1, m_indexStream);

    // store offset of beginning of data
    m_dataBeginOffset = ftell64(m_indexStream);

    // return success/failure of write
    return ( elementsWritten == 6 );
}

// write index data for all references to new index file
bool BamToolsIndex::WriteAllReferences(void) {

    size_t elementsWritten = 0;

    // write number of references
    int32_t numReferences = (int32_t)m_indexData.size();
    if ( m_isBigEndian ) SwapEndian_32(numReferences);
    elementsWritten += fwrite(&numReferences, sizeof(numReferences), 1, m_indexStream);

    // iterate through references in index
    bool refOk = true;
    BamToolsIndexData::const_iterator refIter = m_indexData.begin();
    BamToolsIndexData::const_iterator refEnd  = m_indexData.end();
    for ( ; refIter != refEnd; ++refIter )
	refOk &= WriteReferenceEntry( (*refIter).second );

    return ( (elementsWritten == 1) && refOk );
}

// write current reference index data to new index file
bool BamToolsIndex::WriteReferenceEntry(const BamToolsReferenceEntry& refEntry) {

    size_t elementsWritten = 0;

    // write number of offsets listed for this reference
    uint32_t numOffsets = refEntry.Offsets.size();
    if ( m_isBigEndian ) SwapEndian_32(numOffsets);
    elementsWritten += fwrite(&numOffsets, sizeof(numOffsets), 1, m_indexStream);

    // iterate over offset entries
    bool entriesOk = true;
    vector<BamToolsIndexEntry>::const_iterator offsetIter = refEntry.Offsets.begin();
    vector<BamToolsIndexEntry>::const_iterator offsetEnd  = refEntry.Offsets.end();
    for ( ; offsetIter != offsetEnd; ++offsetIter )
	entriesOk &= WriteIndexEntry( (*offsetIter) );

    return ( (elementsWritten == 1) && entriesOk );
}

// write current index offset entry to new index file
bool BamToolsIndex::WriteIndexEntry(const BamToolsIndexEntry& entry) {

    // copy entry data
    int32_t maxEndPosition = entry.MaxEndPosition;
    int64_t startOffset    = entry.StartOffset;
    int32_t startPosition  = entry.StartPosition;

    // swap endian-ness if necessary
    if ( m_isBigEndian ) {
	SwapEndian_32(maxEndPosition);
	SwapEndian_64(startOffset);
	SwapEndian_32(startPosition);
    }

    // write the reference index entry
    size_t elementsWritten = 0;
    elementsWritten += fwrite(&maxEndPosition, sizeof(maxEndPosition), 1, m_indexStream);
    elementsWritten += fwrite(&startOffset,    sizeof(startOffset),    1, m_indexStream);
    elementsWritten += fwrite(&startPosition,  sizeof(startPosition),  1, m_indexStream);
    return ( elementsWritten == 3 );
}
