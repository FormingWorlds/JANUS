You are reviewing the JANUS atmosphere model's documentation against its source code, looking for places where they disagree.

JANUS computes a 1D temperature-pressure profile for a planetary atmosphere. The documentation below describes the physical model, including numbered equations that in several places cite specific equations from published papers (e.g. "Graham et al. 2021 Eq. 15"). The source code below implements this model.

Your job: compare the two and report every place they disagree, plus every place the documentation is silent about something the code does. Respond with a single JSON object matching the required schema, containing everything you find (an empty `findings` array is a valid and expected result if you find nothing).

For each finding:
- Quote the *exact* doc passage in `doc_excerpt` and the *exact* code in `code_excerpt` — do not paraphrase either side. This lets a human reviewer verify your claim without re-deriving it themselves.
- Use `code_lines` referring to the line numbers shown in the numbered source listing below (e.g. "L42-L58"). Those `N: ` prefixes are for reference only: do not include them in `code_excerpt`. The documentation is shown unnumbered, exactly as it is in the file.
- For a `gap`, set `doc_excerpt` to the closest related doc passage if there is one, or to an empty string if there is none.
- Classify as `type: "inconsistency"` when the doc and code both address the same thing but disagree (a different equation, a different sign, a different coefficient, a different assumption). Set `severity: "serious"` if the disagreement would change model output or physical interpretation, `"minor"` if it's a simplification, notation difference, or something that doesn't change results.
- Classify as `type: "gap"` when the code does something the docs don't mention at all. Do not set `severity` for gaps.
- Only include `suggested_fix` when you can give a safe, literal, verbatim text replacement within `doc_file` — `old_text` must be copied character-for-character from your `doc_excerpt`, must occur exactly once in the doc file (include enough surrounding words to make it unique), and must differ from `new_text`, its replacement. If the right fix isn't a simple text swap (it needs new equations, new sections, or a judgement call about which reference is authoritative), omit `suggested_fix` entirely rather than guessing. An incorrect physics claim is worse than no suggested fix.
- Never invent a citation, equation, or numerical value that doesn't appear in the text below. If you're not sure whether something is actually a disagreement (e.g. you can't tell if a variable name means the same thing in both places), say so in `description` and still report it, but don't assert certainty you don't have.

--- DOCUMENTATION ---
{{DOCS}}

--- SOURCE CODE ---
{{SOURCE}}
