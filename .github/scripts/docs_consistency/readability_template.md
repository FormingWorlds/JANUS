You are editing proposed fixes to the documentation of the JANUS atmosphere model. Each fix replaces `old_text` in the doc with `new_text`. The `new_text` is already correct; your only job is to make it read clearly for a scientist who knows atmospheric physics but has not read the JANUS source code.

For each fix, return its `id` and a `new_text`. If the given `new_text` already reads clearly, return it unchanged.

What to improve:
- One idea per sentence. Split a sentence that strings several clauses together with commas, or that holds a long aside in brackets, into shorter sentences.
- Keep a list or an "either ... or" structure intact: if a long clause breaks it, move that clause into its own sentence.
- Plain words over code-speak: describe what the model does ("the tropopause search is skipped"), not the code's control flow ("the step is skipped", "the branch is taken").
- Name code identifiers (in backticks) only for settings a user can change or for functions a reader is pointed to; never use one in place of a physical description.
- Make directions explicit in a pressure column: say "at higher pressure" or "deeper" rather than a bare "below" or "above" when either reading is possible.
- Do not restate what `old_text` already said elsewhere in the paragraph (see `doc_excerpt`).

What you must keep:
- Every statement, condition, number, unit, symbol, code identifier and citation in the given `new_text`. Add no new claims, and drop none. The `code_excerpt` and `description` show what the text must stay true to; do not use them to add anything.
- The doc's notation: every equation, symbol or variable as LaTeX inside `$...$`, written exactly as in the given `new_text`.
- The fit with its surroundings: `new_text` replaces `old_text` inside `doc_excerpt`, so it must start and end so that the text around it still reads as one passage.
- Line breaks and indentation: where `old_text` is wrapped by hand (e.g. the indented text of a `!!! note`), wrap `new_text` at the same width and indentation; it may end up a line longer or shorter. Keep every paragraph break, list item and other line-start marker.
- No em-dashes or en-dashes; use a colon, semicolon, comma or parentheses instead.

The fixes, as JSON:

{{FIXES}}
