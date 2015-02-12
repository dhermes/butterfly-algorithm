import _ast
import ast
import os
import subprocess


ACCEPTED_NODE_TYPES = (_ast.Import, _ast.ImportFrom, _ast.FunctionDef)
GENERATED_FILENAME = 'butterfly.py'
DESIRED_FILENAME = 'butterfly_algorithm.py'
UNUSED_MODULES = ('make_partition_plots', 'IPython', 'matplotlib',
                  'time', 'os', 'sympy')


def generate_file():
    # NOTE: We could probably do this with `import IPython` ... etc.
    #       but the overhead is not likely worth it.
    subprocess.check_output(['ipython', 'nbconvert', '--to',
                             'python', 'butterfly.ipynb'])


def get_tree():
    with open(GENERATED_FILENAME, 'rU') as fh:
        generated_butterfly = fh.read()

    file_lines = generated_butterfly.split('\n')
    tree = ast.parse(generated_butterfly)
    return tree, file_lines


def get_captured_lines(tree):
    captured_imports = []
    captured_functions = []
    for i, node in enumerate(tree.body):
        if isinstance(node, ACCEPTED_NODE_TYPES):
            # tree.body should only be top level.
            if node.col_offset != 0:
                raise ValueError('Node is not top-level', node)
            # NOTE: This may fail if `node` is the last entry.
            next_node = tree.body[i + 1]

            section = (node.lineno - 1, next_node.lineno - 1)
            if isinstance(node, _ast.FunctionDef):
                # The `get_time` function requires globals from the
                # notebook which we don't have and we don't need to use
                # the `custom_update` function.
                if node.name not in ('get_time', 'custom_update'):
                    captured_functions.append(section)
            elif isinstance(node, _ast.ImportFrom):
                if node.module not in UNUSED_MODULES:
                    captured_imports.append(section)
            else:
                if (len(node.names) == 1 and
                    node.names[0].name not in UNUSED_MODULES):
                    captured_imports.append(section)

    return captured_imports, captured_functions


def write_import(lines, fh):
    for line in lines:
        # Don't include top level comments.
        if line.startswith('#'):
            continue
        # Don't write blank lines.
        if not line.strip():
            continue

        # If the line is accepted write to the file.
        fh.write(line + '\n')


def write_function(lines, fh):
    true_last_line = len(lines) - 1
    for last_line in xrange(len(lines) - 1, -1, -1):
        if lines[last_line].startswith('#'):
            continue
        elif not lines[last_line].strip():
            continue
        else:
            true_last_line = last_line
            break

    for line in lines[:true_last_line + 1]:
        fh.write(line + '\n')


def rewrite_file():
    tree, file_lines = get_tree()
    captured_imports, captured_functions = get_captured_lines(tree)

    with open(DESIRED_FILENAME, 'w') as fh:
        for begin, end in captured_imports:
            write_import(file_lines[begin:end], fh)

        for begin, end in captured_functions:
            # Two newlines between every function
            fh.write('\n\n')
            write_function(file_lines[begin:end], fh)


def main():
    generate_file()
    print 'Writing functions to: %r' % (DESIRED_FILENAME,)
    rewrite_file()
    print 'Removing: %r' % (GENERATED_FILENAME,)
    os.remove(GENERATED_FILENAME)


if __name__ == '__main__':
    main()
