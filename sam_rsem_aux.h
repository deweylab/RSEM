#ifndef SAM_RSEM_AUX_H_
#define SAM_RSEM_AUX_H_

#include<cstdio>
#include<cstring>
#include<stdint.h>

#include "sam/bam.h"

// dwt: duplicate without text
bam_header_t *bam_header_dwt(const bam_header_t *ori_h)
{
	bam_header_t *h;

	h = bam_header_init();
	h->n_targets = ori_h->n_targets;
	h->target_len = (uint32_t*)calloc(h->n_targets, 4);
	h->target_name = (char**)calloc(h->n_targets, sizeof(char*));
	for (int i = 0; i < h->n_targets; i++) {
		h->target_len[i] = ori_h->target_len[i];
		h->target_name[i] = strdup(ori_h->target_name[i]);
	}

	return h;
}

void append_header_text(bam_header_t *header, const char* text, int len)
{
	int x = header->l_text + 1;
	int y = header->l_text + len + 1; // 1 byte null
	if (text == 0) return;
	kroundup32(x);
	kroundup32(y);
	if (x < y) header->text = (char*)realloc(header->text, y);
	strncpy(header->text + header->l_text, text, len); // we cannot use strcpy() here.
	header->l_text += len;
	header->text[header->l_text] = 0;
}

void expand_data_size(bam1_t *b) {
	if (b->m_data < b->data_len) {
		b->m_data = b->data_len;
		kroundup32(b->m_data);
		b->data = (uint8_t*)realloc(b->data, b->m_data);
	}
}

#endif /* SAM_RSEM_AUX_H_ */
