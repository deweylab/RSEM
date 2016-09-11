// ***************************************************************************
// BamWriter.cpp (c) 2009 Michael Str�mberg, Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 22 November 2010 (DB)
// ---------------------------------------------------------------------------
// Provides the basic functionality for producing BAM files
// ***************************************************************************

#include <BamWriter.h>
#include <BamWriter_p.h>
using namespace BamTools;
using namespace BamTools::Internal;

#include <iostream>
using namespace std;

// constructor
BamWriter::BamWriter(void) {
    d = new BamWriterPrivate;
}

// destructor
BamWriter::~BamWriter(void) {
    delete d;
    d = 0;
}

// closes the alignment archive
void BamWriter::Close(void) {
    d->Close();
}

// opens the alignment archive
bool BamWriter::Open(const string& filename,
                     const string& samHeader,
                     const RefVector& referenceSequences,
                     bool isWriteUncompressed)
{
    return d->Open(filename, samHeader, referenceSequences, isWriteUncompressed);
}

// saves the alignment to the alignment archive
void BamWriter::SaveAlignment(const BamAlignment& al) {
    d->SaveAlignment(al);
}
