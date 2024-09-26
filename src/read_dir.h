#ifndef READ_DIR_H
#define READ_DIR_H
#include "integers.h"

typedef struct FileInfo {
    char **dirnames, **filenames;
    int32 nfiles;
    int32 ndirs;
} FileInfo;

#define MAXPATHLEN 512
extern char cur_dir[MAXPATHLEN];
extern FileInfo my_ff;

void read_dir_free_finfo(FileInfo *ff);
int32 read_dir_get_fileinfo_tab(char *wild, char *direct, FileInfo *ff, char *wild2);
int32 read_dir_get_fileinfo(char *wild, char *direct, FileInfo *ff);
int32 read_dir_fil_count(char *direct, int32 *ndir, int32 *nfil, char *wild, int32 *mld,
                int32 *mlf);
int32 read_dir_change_dir(char *path);
int32 read_dir_get_directory(char *direct);
int32 read_dir_is_directory(char *root, char *path);
int32 read_dir_wild_match(char *string, char *pattern);

#endif
