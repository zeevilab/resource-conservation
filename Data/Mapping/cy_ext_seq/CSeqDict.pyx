from libc.stdio cimport *
from libcpp.string cimport string
from libcpp.unordered_map cimport unordered_map
from libcpp.pair cimport pair
from libc.stdlib cimport malloc, free
cimport cython
from cpython.mem cimport PyMem_Malloc, PyMem_Free
import ujson
from os.path import splitext

cdef char *basemap = <char *>PyMem_Malloc(255 * cython.sizeof(char))
cdef int c_i = 0
while c_i < 255:
    basemap[c_i] = 'N'
    c_i += 1
basemap[65] = 'T'
basemap[67] = 'G'
basemap[71] = 'C'
basemap[84] = 'A'
basemap[97] = 'T'
basemap[99] = 'G'
basemap[103] = 'C'
basemap[116] = 'A'
basemap[78] = 'N'

PLUS = '+'
# cdef class SeqIDLL:
#     cdef char *name
#     cdef SeqIDLL newer
#     cdef SeqIDLL older

# cpdef cppclass SeqIDLL:
#     string name
#     SeqIDLL *newer
#     SeqIDLL *older

cdef SeqIDLL * newSeqID(string seqid):
    cdef SeqIDLL *new_ = new SeqIDLL();
    new_.name = seqid
    new_.newer = EMPTY
    new_.older = EMPTY
    return new_
    
cdef struct TestS:
    string name
 
cpdef int test(string n) except -1:
    cdef TestS cur
    cur.name = n
    return 0
    
  
# ctypedef s_SeqIDLL SeqIDLL

cdef char * revcomp(unsigned char *seq, long seqlen):
    cdef char *ret = <char *> PyMem_Malloc(seqlen * cython.sizeof(char) + 1) 
    if ret == NULL:
        print 'mem in revcomp'
    cdef long i = 0
    while i < seqlen:
        ret[i] = basemap[seq[seqlen-i-1]]
        i+=1
    ret[seqlen] = '\0'
    return ret

cdef SeqIDLL *EMPTY

cdef class CSeqDict:
#     cdef unordered_map[string, char *] _dct
#     cdef unordered_map[string, char *] _revdct
#     cdef unordered_map[string, long] _lendct
#     cdef unordered_map[string, long] _idx 
#     cdef unordered_map[string, SeqIDLL *] _cppcq
#     cdef long _cursum, _lim
#     cdef FILE* _file
#     cdef SeqIDLL *newest
#     cdef SeqIDLL *oldest
    
    def __init__(self, pth, lim = None):
        with open(splitext(pth)[0] + '.lens') as ol:
            for k, v in ujson.load(ol).items():
                self._lendct[k] = v 
        with open(splitext(pth)[0] + '.idx') as o:
            for k, v in ujson.load(o).items():
                self._idx[k] = v
        self._cursum = 0
#         self._cacheq = {}
        self._lim = lim
        self._file = fopen(pth, "rb")
        if self._file == NULL:
            raise Exception(2, "No such file or directory: '%s'" % pth)
    
    def get_lendct(self):
        return self._lendct
    
    def __dealloc__(self):
        fclose(self._file)
    
    def _get_rev(self, seqid):
        return self.get_rev(seqid)
    
    def _get_seq(self, seqid):
        return self.get_seq(seqid, 1)
    
    def get_sequence(self, gid, strand, pos, lngt):
        ln = self.get_len(gid)
        if strand == PLUS:
            sq, st, en, strnd =  self.get_seq(gid, 1), pos-1, pos+lngt-1, 1 
        else:
            sq = self.get_rev(gid)
            st, en, strnd =  ln-pos-lngt+1, ln-pos+1, -1
        return sq, st%ln, en%ln, strnd
    
    cdef char * get_seq(self, string seqid, int callcache) except NULL:
        cdef char *sq
        if callcache:
            self._add_to_cacheq(seqid)
        if self._dct.count(seqid) > 0:
            return self._dct[seqid]
        else:
            sq = self._load_seq(seqid)
            self._dct[seqid] = sq
#             PyMem_Free(sq)
            self._cursum += self._lendct[seqid]
            self._clean_mem()
            return self._dct[seqid]
    
    cdef char * get_rev(self, string seqid) except NULL:
        cdef long slen
        cdef char *sq
        self._add_to_cacheq(seqid)
        if self._revdct.count(seqid) > 0:
            return self._revdct[seqid]
        else:
            slen = self._lendct[seqid]
            sq = revcomp(<unsigned char *>self.get_seq(seqid, 0), slen)
            self._revdct[seqid] = sq
#             PyMem_Free(sq)
            self._cursum += slen
            self._clean_mem()
            return self._revdct[seqid]
        
#     cpdef int get_len(self, char *seqid):
    cpdef long get_len(self, string seqid):
        return self._lendct[seqid]
    
    cpdef long get_cursum(self):
        return self._cursum

    
        

    cdef int _add_to_cacheq(self, string seqid) except -1:
        cdef SeqIDLL *cur
        if self._cppcq.count(seqid) > 0:
            cur = self._cppcq[seqid]
            if cur.newer == EMPTY:
                return 0
            if cur.older == EMPTY:
                cur.newer.older = EMPTY
                self.oldest = cur.newer
            else:
                cur.older.newer = cur.newer
                cur.newer.older = cur.older
            self.newest.newer = cur
            cur.newer = EMPTY
            cur.older = self.newest
            self.newest = cur
        else:
            cur = newSeqID(seqid)
            cur.newer = EMPTY
            cur.older = EMPTY
            if self.newest == NULL:
                self.newest = cur
                self.oldest = cur
            else:
                self.newest.newer = cur
                cur.newer = EMPTY
                cur.older = self.newest
                self.newest = cur
            self._cppcq[seqid] = cur
        return 0
    
    cdef char * _load_seq(self, string seqid):
        cdef long pos = self._idx[seqid]
        cdef char * line = NULL
        cdef size_t l = 0
        fseek(self._file, pos, SEEK_SET)
        getline(&line, &l, self._file)
         
        cdef long lseq = self._lendct[seqid]
        cdef long ind = lseq - 2
        while line[ind] != 10 and line[ind] != 0:
            ind += 1
        line[ind] = '\0' 
        return line 
    
    cdef int _clean_mem(self) except -1:
        cdef string seqtodl
        cdef SeqIDLL *new_oldest = NULL
        while self._cursum > self._lim:
            seqtodl = self.oldest.name
            self._cppcq.erase(seqtodl)
            PyMem_Free(self._dct[seqtodl])
            self._dct.erase(seqtodl)
            self._cursum -= self._lendct[seqtodl]
            if self._revdct.count(seqtodl) > 0:
                PyMem_Free(self._revdct[seqtodl])
                self._revdct.erase(seqtodl)
                self._cursum -= self._lendct[seqtodl]
            if self.oldest.newer != NULL:
                new_oldest = self.oldest.newer
                new_oldest.older = EMPTY
            PyMem_Free(self.oldest)
            self.oldest = new_oldest
            if new_oldest == NULL:
                self.newest = NULL
        return 0
#     cdef int _add_to_cacheq(self, char *seqid) except -1:
#         self._cacheq[seqid] = self._count
#         self._count += 1
#         
#     def _clean_mem(self):
#         while self._cursum > self._lim:
#             seqtodl = min(self._cacheq, key = self._cacheq.get)
#             del self._cacheq[seqtodl]
#             del self._dct[seqtodl]
#             self._cursum -= self._lendct[seqtodl] * 2
#             try: 
#                 del self._revdct[seqtodl]
#                 self._cursum -= self._lendct[seqtodl] * 2
#             except: pass
#         
#     cdef int _add_to_cacheq(self, char *seqid) except -1:
#         cdef SeqIDLL cur
#         if seqid in self._cacheq:
#             cur = self._cacheq[seqid]
#             if cur.newer == EMPTY:
#                 return 0
#             if cur.older == EMPTY:
#                 cur.newer.older = EMPTY
#                 self.oldest = cur.newer
#             else:
#                 cur.older.newer = cur.newer
#                 cur.newer.older = cur.older
#             self.newest.newer = cur
#             cur.newer = EMPTY
#             cur.older = self.newest
#             self.newest = cur
#         else:
#             cur = SeqIDLL()
#             cur.name = seqid
#             cur.newer = EMPTY
#             cur.older = EMPTY
#             if self.newest is None:
#                 self.newest = cur
#                 self.oldest = cur
#             else:
#                 self.newest.newer = cur
#                 cur.newer = EMPTY
#                 cur.older = self.newest
#                 self.newest = cur
#             self._cacheq[seqid] = cur
#         return 0
#     
#     cdef char * _load_seq(self, char *seqid):
#         cdef long pos = self._idx[seqid]
#         cdef char * line = NULL
#         cdef size_t l = 0
#         fseek(self._file, pos, SEEK_SET)
#         getline(&line, &l, self._file)
#         if l > self._lendct[seqid] + 1 and line[l-2] == 10:
#             line[l-2] = '\0' 
#         return line 
#     
#     cdef int _clean_mem(self) except -1:
#         cdef char *seqtodl
#         cdef SeqIDLL new_oldest
#         while self._cursum > self._lim:
#             seqtodl = self.oldest.name
#             del self._cacheq[seqtodl]
#             del self._dct[seqtodl]
#             self._cursum -= self._lendct[seqtodl]
#             if seqtodl in self._revdct:
#                 del self._revdct[seqtodl]
#                 self._cursum -= self._lendct[seqtodl]
#             new_oldest = self.oldest.newer
#             new_oldest.older = EMPTY
#             self.oldest = new_oldest
#         return 0
            
# from libc.stdio cimport *
# from libcpp.string cimport string
# from libcpp.unordered_map cimport unordered_map
# from libcpp.unordered_set cimport unordered_set 
# from libcpp.pair cimport pair
# 
# cimport cython
# from cpython.mem cimport PyMem_Malloc, PyMem_Free
# import ujson
# from os.path import splitext
# 
# cdef char *basemap = <char *>PyMem_Malloc(255 * cython.sizeof(char))
# cdef int c_i = 0
# while c_i < 255:
#     basemap[c_i] = '\0'
#     c_i += 1
# basemap[65] = 'T'
# basemap[67] = 'G'
# basemap[71] = 'C'
# basemap[84] = 'A'
# basemap[97] = 'T'
# basemap[99] = 'G'
# basemap[103] = 'C'
# basemap[116] = 'A'
# basemap[78] = 'N'
# 
# # cdef class SeqIDLL:
# #     cdef char *name
# #     cdef SeqIDLL newer
# #     cdef SeqIDLL older
# 
# cdef struct SeqIDLL:
#     string *name
#     SeqIDLL *newer
#     SeqIDLL *older
# 
# 
# # cdef class TestS:
# #     cdef char *name
# # 
# # cdef dict test_dict = {}
# # cpdef int test(char* n) except -1:
# #     cdef TestS cur = TestS()
# #     cur.name = n
# #     print cur.name
# #     test_dict[n] = cur
# #     print cur.name
# #     return 0
#     
#   
# # ctypedef s_SeqIDLL SeqIDLL
# 
# cdef char * revcomp(unsigned char *seq, int seqlen):
#     cdef char *ret = <char *> PyMem_Malloc(seqlen * cython.sizeof(char)) 
#     cdef int i = 0
#     while i < seqlen:
#         ret[i] = basemap[seq[seqlen-i-1]]
#         i+=1
#     return ret
# 
# cdef SeqIDLL *EMPTY
# 
# cdef class CSeqDict:
#     cdef dict _dct, _revdct, _lendct, _idx#, _cacheq
#     cdef unordered_map[string, SeqIDLL *] _cppcq
#     cdef unordered_set[string] _revs
#     cdef long _cursum, _lim
#     cdef FILE* _file
#     cdef SeqIDLL *newest
#     cdef SeqIDLL *oldest
#     
#     def __init__(self, pth, lim = None):
#         self._dct = {}
#         self._revdct = {}
#         with open(splitext(pth)[0] + '.lens') as ol:
#             self._lendct = ujson.load(ol)
#         with open(splitext(pth)[0] + '.idx') as o:
#             self._idx = ujson.load(o)
#         self._cursum = 0
# #         self._cacheq = {}
#         self._lim = lim
#         self._file = fopen(pth, "rb")
#         if self._file == NULL:
#             raise Exception(2, "No such file or directory: '%s'" % pth)
#     
#     def __dealloc__(self):
#         fclose(self._file)
#     
# #     cpdef char * get_seq(self, char *seqid, int callcache = 1):
#     cpdef char * get_seq(self, string seqid, int callcache = 1):
#         cdef char *sq
# #         if callcache:
# #             self._add_to_cacheq(seqid)
# #         if seqid in self._dct:
#         if self._cppcq.count(seqid) > 0:
#             if callcache:
#                 self._add_to_cacheq(seqid)
#             return self._dct[seqid]
#         else:
#             if callcache:
#                 self._add_to_cacheq(seqid)
#             sq = self._load_seq(seqid)
#             self._dct[seqid] = sq
#             PyMem_Free(sq)
#             self._cursum += self._lendct[seqid]
#             self._clean_mem()
#             return self._dct[seqid]
#     
# #     cpdef char * get_rev(self, char *seqid):    
#     cpdef char * get_rev(self, string seqid):
#         cdef int slen
#         cdef char *sq
#         self._add_to_cacheq(seqid)
#         if self._revs.count(seqid) > 0:
# #         if seqid in self._revdct:
#             return self._revdct[seqid]
#         else:
#             slen = self._lendct[seqid]
#             sq = revcomp(<unsigned char *>self.get_seq(seqid, 0), slen) 
#             self._revdct[seqid] = sq
#             PyMem_Free(sq)
#             self._revs.insert(seqid)
#             self._cursum += slen
#             self._clean_mem()
#             return self._revdct[seqid]
#         
# #     cpdef int get_len(self, char *seqid):
#     cpdef int get_len(self, string seqid):
#         return self._lendct[seqid]
#     
#     cpdef int get_cursum(self):
#         return self._cursum
#     
#     cdef int _add_to_cacheq(self, string seqid) except -1:
#         cdef SeqIDLL *cur
#         if self._cppcq.count(seqid) > 0:
#             cur = self._cppcq[seqid]
#             if cur.newer == EMPTY:
#                 return 0
#             if cur.older == EMPTY:
#                 cur.newer.older = EMPTY
#                 self.oldest = cur.newer
#             else:
#                 cur.older.newer = cur.newer
#                 cur.newer.older = cur.older
#             self.newest.newer = cur
#             cur.newer = EMPTY
#             cur.older = self.newest
#             self.newest = cur
#         else:
#             cur = <SeqIDLL*>PyMem_Malloc(cython.sizeof(SeqIDLL))
#             cur.name = &seqid
#             cur.newer = EMPTY
#             cur.older = EMPTY
#             if self.newest == NULL:
#                 self.newest = cur
#                 self.oldest = cur
#             else:
#                 self.newest.newer = cur
#                 cur.newer = EMPTY
#                 cur.older = self.newest
#                 self.newest = cur
#             self._cppcq[seqid] = cur
#         return 0
#     
#     cdef char * _load_seq(self, string seqid):
#         cdef long pos = self._idx[seqid]
#         cdef char * line = NULL
#         cdef size_t l = 0
#         fseek(self._file, pos, SEEK_SET)
#         getline(&line, &l, self._file)
#         if l > self._lendct[seqid] + 1 and line[l-2] == 10:
#             line[l-2] = '\0' 
#         return line 
#     
#     cdef int _clean_mem(self) except -1:
#         cdef string seqtodl
#         cdef SeqIDLL *new_oldest
#         while self._cursum > self._lim:
#             seqtodl = self.oldest.name[0]
#             self._cppcq.erase(seqtodl)
#             del self._dct[seqtodl]
#             self._cursum -= self._lendct[seqtodl]
#             if self._revs.count(seqtodl) > 0:
# #             if seqtodl in self._revdct:
#                 del self._revdct[seqtodl]
#                 self._revs.erase(seqtodl)
#                 self._cursum -= self._lendct[seqtodl]
#             new_oldest = self.oldest.newer
#             new_oldest.older = EMPTY
#             PyMem_Free(self.oldest)
#             self.oldest = new_oldest
#         return 0
# #     cdef int _add_to_cacheq(self, char *seqid) except -1:
# #         self._cacheq[seqid] = self._count
# #         self._count += 1
# #         
# #     def _clean_mem(self):
# #         while self._cursum > self._lim:
# #             seqtodl = min(self._cacheq, key = self._cacheq.get)
# #             del self._cacheq[seqtodl]
# #             del self._dct[seqtodl]
# #             self._cursum -= self._lendct[seqtodl] * 2
# #             try: 
# #                 del self._revdct[seqtodl]
# #                 self._cursum -= self._lendct[seqtodl] * 2
# #             except: pass
# #         
# #     cdef int _add_to_cacheq(self, char *seqid) except -1:
# #         cdef SeqIDLL cur
# #         if seqid in self._cacheq:
# #             cur = self._cacheq[seqid]
# #             if cur.newer == EMPTY:
# #                 return 0
# #             if cur.older == EMPTY:
# #                 cur.newer.older = EMPTY
# #                 self.oldest = cur.newer
# #             else:
# #                 cur.older.newer = cur.newer
# #                 cur.newer.older = cur.older
# #             self.newest.newer = cur
# #             cur.newer = EMPTY
# #             cur.older = self.newest
# #             self.newest = cur
# #         else:
# #             cur = SeqIDLL()
# #             cur.name = seqid
# #             cur.newer = EMPTY
# #             cur.older = EMPTY
# #             if self.newest is None:
# #                 self.newest = cur
# #                 self.oldest = cur
# #             else:
# #                 self.newest.newer = cur
# #                 cur.newer = EMPTY
# #                 cur.older = self.newest
# #                 self.newest = cur
# #             self._cacheq[seqid] = cur
# #         return 0
# #     
# #     cdef char * _load_seq(self, char *seqid):
# #         cdef long pos = self._idx[seqid]
# #         cdef char * line = NULL
# #         cdef size_t l = 0
# #         fseek(self._file, pos, SEEK_SET)
# #         getline(&line, &l, self._file)
# #         if l > self._lendct[seqid] + 1 and line[l-2] == 10:
# #             line[l-2] = '\0' 
# #         return line 
# #     
# #     cdef int _clean_mem(self) except -1:
# #         cdef char *seqtodl
# #         cdef SeqIDLL new_oldest
# #         while self._cursum > self._lim:
# #             seqtodl = self.oldest.name
# #             del self._cacheq[seqtodl]
# #             del self._dct[seqtodl]
# #             self._cursum -= self._lendct[seqtodl]
# #             if seqtodl in self._revdct:
# #                 del self._revdct[seqtodl]
# #                 self._cursum -= self._lendct[seqtodl]
# #             new_oldest = self.oldest.newer
# #             new_oldest.older = EMPTY
# #             self.oldest = new_oldest
# #         return 0
            
            
    
    
            
    
    
