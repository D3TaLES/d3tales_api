#!/usr/bin/env bash
if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
    bash "${BASH_SOURCE[0]}" "$@"
    return $?
fi

set -euo pipefail

version="0.2"
author="Duke, Rebekah"
projectname="D<sup>3</sup>TaLES API"

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd -- "$script_dir/.." && pwd)"
build_dir="$(mktemp -d)"
source_dir="$build_dir/source"
html_dir="$build_dir/html"
output_dir="$repo_root/docs"

cleanup() {
    rm -rf -- "$build_dir"
}
trap cleanup EXIT

mkdir -p "$source_dir" "$html_dir" "$output_dir"
output_dir="$(cd -- "$output_dir" && pwd)"

export DB_INFO_FILE="$repo_root/db_info_ex.json"
export D3TALES_API_ROOT="$repo_root"

sphinx-apidoc --force \
    --no-toc \
    --doc-project "$projectname" \
    --doc-author "$author" \
    --doc-version "$version" \
    --output-dir "$source_dir" \
    "$repo_root/d3tales_api" \
    "$repo_root/d3tales_api/database_info" \
    "$repo_root/d3tales_api/Processors/parser_nlp.py"

cp "$script_dir/conf.py" "$script_dir/index.rst" "$script_dir"/*.md "$script_dir"/*.png "$source_dir/"

sphinx-build -W -b html -E -a -d "$build_dir/doctrees" "$source_dir" "$html_dir"
rm -f -- "$html_dir/.buildinfo"
find "$output_dir" -mindepth 1 -maxdepth 1 -exec rm -rf -- {} +
cp -a "$html_dir/." "$output_dir/"
touch "$output_dir/.nojekyll"

echo "Documentation generated in $output_dir"
