// ***************************************************************************
// BamReader.cpp (c) 2009 Derek Barnett, Michael Str�mberg
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 22 November 2010 (DB)
// ---------------------------------------------------------------------------
// Provides the basic functionality for reading BAM files
// ***************************************************************************

#include <BamReader.h>
#include <BamReader_p.h>
using namespace BamTools;
using namespace BamTools::Internal;

#include <algorithm>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>
using namespace std;

// constructor
BamReader::BamReader(void) {
    d = new BamReaderPrivate(this);
}

// destructor
BamReader::~BamReader(void) {
    delete d;
    d = 0;
}

// file operations
void BamReader::Close(void) { d->Close(); }
bool BamReader::HasIndex(void) const { return d->HasIndex; }
bool BamReader::IsIndexLoaded(void) const { return HasIndex(); }
bool BamReader::IsOpen(void) const { return d->mBGZF.IsOpen; }
bool BamReader::Jump(int refID, int position)  { return d->SetRegion( BamRegion(refID, position) ); }
bool BamReader::Open(const std::string& filename,
                     const std::string& indexFilename,
                     const bool lookForIndex,
                     const bool preferStandardIndex)
{
    return d->Open(filename, indexFilename, lookForIndex, preferStandardIndex);
}
bool BamReader::Rewind(void) { return d->Rewind(); }
bool BamReader::SetRegion(const BamRegion& region) { return d->SetRegion(region); }
bool BamReader::SetRegion(const int& leftRefID, const int& leftBound, const int& rightRefID, const int& rightBound) {
    return d->SetRegion( BamRegion(leftRefID, leftBound, rightRefID, rightBound) );
}

// access alignment data
bool BamReader::GetNextAlignment(BamAlignment& bAlignment) { return d->GetNextAlignment(bAlignment); }
bool BamReader::GetNextAlignmentCore(BamAlignment& bAlignment) { return d->GetNextAlignmentCore(bAlignment); }

// access auxiliary data
const string BamReader::GetHeaderText(void) const { return d->GetHeaderText(); }
int BamReader::GetReferenceCount(void) const { return d->References.size(); }
const RefVector& BamReader::GetReferenceData(void) const { return d->References; }
int BamReader::GetReferenceID(const string& refName) const { return d->GetReferenceID(refName); }
const std::string BamReader::GetFilename(void) const { return d->Filename; }

// index operations
bool BamReader::CreateIndex(bool useStandardIndex) { return d->CreateIndex(useStandardIndex); }
void BamReader::SetIndexCacheMode(const BamIndex::BamIndexCacheMode mode) { d->SetIndexCacheMode(mode); }
