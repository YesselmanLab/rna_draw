/*
 * vrna_compat.c -- minimal ViennaRNA compatibility shim for the vendored
 * RNApuzzler / RNAturtle layout core.
 *
 * PURPOSE
 * -------
 * The vendored layout objects (vendor/RNApuzzler/RNApuzzler.c and
 * RNAturtle.c) reference EXACTLY TWO ViennaRNA library symbols that live in
 * libRNA.a:
 *
 *   - vrna_alloc   (upstream src/ViennaRNA/utils/utils.c)  -- a calloc wrapper
 *   - vrna_ptable  (upstream src/ViennaRNA/structures/structure_pairtable.c),
 *                  which reduces to vrna_ptable_from_string + a file-static
 *                  extract_pairs helper.
 *
 * Vendoring these ~self-contained definitions lets `_vienna_layout` compile
 * and link the puzzler + turtle cores WITHOUT linking libRNA.a at all --
 * i.e. rna_draw has NO ViennaRNA runtime dependency for layout. The output
 * is byte-identical to linking libRNA.a for these two symbols (verified:
 * these functions are pure, deterministic, and copied faithfully from the
 * upstream 2.7.0 sources).
 *
 * The only other libRNA.a-reachable symbol is the log/error reporter that
 * vrna_alloc and vrna_ptable_from_string call ONLY on allocation failure or
 * malformed/oversized input. The _vienna_layout binding pre-validates every
 * structure (well-nestedness, non-empty, no empty loop) before calling in,
 * so those branches are unreachable in practice; we vendor `vrna_log` as a
 * trivial stderr stub to satisfy the linker.
 *
 * PROVENANCE
 * ----------
 * All definitions below are copied faithfully from ViennaRNA 2.7.0
 * (https://github.com/ViennaRNA/ViennaRNA/releases/tag/v2.7.0):
 *   - vrna_alloc                  <- src/ViennaRNA/utils/utils.c
 *   - vrna_ptable / from_string   <- src/ViennaRNA/structures/structure_pairtable.c
 *   - extract_pairs               <- src/ViennaRNA/structures/structure_pairtable.c
 * The upstream PUBLIC/PRIVATE/INLINE decorations are collapsed to plain C
 * linkage / `static`; the `vrna_log_*` printf-macros are routed through the
 * `vrna_log` stub. License: upstream ViennaRNA is GPLv2+; used here with the
 * project owner's permission (see vendor/README.md).
 */

#include <errno.h>
#include <limits.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Declarations we are defining (keeps signatures in lockstep with the
 * ViennaRNA headers the vendored C and the binding compile against).
 * ViennaRNA/structures/dotbracket.h defines the VRNA_BRACKETS_* bitmask but
 * transitively drags in the fold-compound / probability headers (FLT_OR_DBL,
 * vrna_fold_compound_t, ...) which are unrelated to layout and don't compile
 * standalone here, so the four bitmask values are reproduced locally instead
 * (verbatim from dotbracket.h, ViennaRNA 2.7.0). */
#include <ViennaRNA/structures/pairtable.h> /* vrna_ptable[_from_string] */
#include <ViennaRNA/utils/log.h>            /* vrna_log_levels_e, vrna_log */

#define VRNA_BRACKETS_ALPHA 4U
#define VRNA_BRACKETS_RND   8U
#define VRNA_BRACKETS_CLY   16U
#define VRNA_BRACKETS_ANG   32U
#define VRNA_BRACKETS_SQR   64U

/* -------------------------------------------------------------------------
 * Message reporter stub.
 *
 * Upstream, the vrna_log_warning()/vrna_log_error() macros funnel into
 * vrna_log(). Those call sites are reached only on OOM or malformed/oversized
 * input, which the _vienna_layout binding rejects up front. A plain stderr
 * reporter is sufficient and keeps the shim free of the ViennaRNA log
 * subsystem (levels, thresholds, callbacks).
 * ------------------------------------------------------------------------- */
void
vrna_log(vrna_log_levels_e  level,
         const char         *file_name,
         int                line_number,
         const char         *format_string,
         ...)
{
  va_list args;

  (void)level;
  fprintf(stderr, "ViennaRNA layout [%s:%d]: ", file_name, line_number);
  va_start(args, format_string);
  vfprintf(stderr, format_string, args);
  va_end(args);
  fputc('\n', stderr);
}


/* -------------------------------------------------------------------------
 * vrna_alloc -- src/ViennaRNA/utils/utils.c (ViennaRNA 2.7.0), faithful copy.
 * ------------------------------------------------------------------------- */
void *
vrna_alloc(size_t size)
{
  void *pointer;

  if ((pointer = (void *)calloc(1, size)) == NULL) {
#ifdef EINVAL
    if (errno == EINVAL) {
      fprintf(stderr, "vrna_alloc: requested size: %lu\n", (unsigned long)size);
      vrna_log_error("Memory allocation failure -> EINVAL");
    }

    if (errno == ENOMEM)
#endif
    vrna_log_error("Memory allocation failure -> no memory");
  }

  return pointer;
}


/* -------------------------------------------------------------------------
 * vrna_ptable and helpers -- src/ViennaRNA/structures/structure_pairtable.c
 * (ViennaRNA 2.7.0), faithful copy. Only the round-bracket path used by
 * vrna_ptable() is exercised by the vendored layout code; the other bracket
 * branches are retained verbatim for fidelity.
 * ------------------------------------------------------------------------- */

/* requires that pt[0] already contains the length of the string! */
static int
extract_pairs(short       *pt,
              const char  *structure,
              const char  *pair)
{
  const char    *ptr;
  char          open, close;
  short         *stack;
  unsigned int  i, j, n;
  int           hx;

  n     = (unsigned int)pt[0];
  stack = (short *)vrna_alloc(sizeof(short) * (n + 1));

  open  = pair[0];
  close = pair[1];

  for (hx = 0, i = 1, ptr = structure; (i <= n) && (*ptr != '\0'); ptr++, i++) {
    if (*ptr == open) {
      stack[hx++] = i;
    } else if (*ptr == close) {
      j = stack[--hx];

      if (hx < 0) {
        vrna_log_warning("%s\nunbalanced brackets '%2s' found while extracting base pairs",
                         structure,
                         pair);
        free(stack);
        return 0;
      }

      pt[i] = j;
      pt[j] = i;
    }
  }

  free(stack);

  if (hx != 0) {
    vrna_log_warning("%s\nunbalanced brackets '%2s' found while extracting base pairs",
                     structure,
                     pair);
    return 0;
  }

  return 1; /* success */
}


short *
vrna_ptable_from_string(const char    *string,
                        unsigned int  options)
{
  char          pairs[3];
  short         *pt;
  unsigned int  i, n;

  n = strlen(string);

  if (n > SHRT_MAX) {
    vrna_log_warning("vrna_ptable_from_string: "
                     "Structure too long to be converted to pair table (n=%d, max=%d)",
                     n,
                     SHRT_MAX);
    return NULL;
  }

  pt    = (short *)vrna_alloc(sizeof(short) * (n + 2));
  pt[0] = (short)n;


  if ((options & VRNA_BRACKETS_RND) &&
      (!extract_pairs(pt, string, "()"))) {
    free(pt);
    return NULL;
  }

  if ((options & VRNA_BRACKETS_ANG) &&
      (!extract_pairs(pt, string, "<>"))) {
    free(pt);
    return NULL;
  }

  if ((options & VRNA_BRACKETS_CLY) &&
      (!extract_pairs(pt, string, "{}"))) {
    free(pt);
    return NULL;
  }

  if ((options & VRNA_BRACKETS_SQR) &&
      (!extract_pairs(pt, string, "[]"))) {
    free(pt);
    return NULL;
  }

  if (options & VRNA_BRACKETS_ALPHA) {
    for (i = 65; i < 91; i++) {
      pairs[0]  = (char)i;
      pairs[1]  = (char)(i + 32);
      pairs[2]  = '\0';
      if (!extract_pairs(pt, string, pairs)) {
        free(pt);
        return NULL;
      }
    }
  }

  return pt;
}


short *
vrna_ptable(const char *structure)
{
  return vrna_ptable_from_string(structure, VRNA_BRACKETS_RND);
}
