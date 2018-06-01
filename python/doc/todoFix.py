## This file is part of OSPREY 3.0
## 
## OSPREY Protein Redesign Software Version 3.0
## Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
## 
## OSPREY is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License version 2
## as published by the Free Software Foundation.
## 
## You should have received a copy of the GNU General Public License
## along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
## 
## OSPREY relies on grants for its development, and since visibility
## in the scientific literature is essential for our success, we
## ask that users of OSPREY cite our papers. See the CITING_OSPREY
## document in this distribution for more information.
## 
## Contact Info:
##    Bruce Donald
##    Duke University
##    Department of Computer Science
##    Levine Science Research Center (LSRC)
##    Durham
##    NC 27708-0129
##    USA
##    e-mail: www.cs.duke.edu/brd/
## 
## <signature of Bruce Donald>, Mar 1, 2018
## Bruce Donald, Professor of Computer Science

import sphinx.ext.todo as todo

# copied from:
# https://raw.githubusercontent.com/sphinx-doc/sphinx/master/sphinx/ext/todo.py

def setup(app):
    
    # override the todo extension function so we can fix it for autodoc TODOs
    # there is apparently some kind of bug with refids
    todo.process_todo_nodes.__code__ = process_todo_nodes.__code__


def process_todo_nodes(app, doctree, fromdocname):
    # type: (Sphinx, nodes.Node, unicode) -> None
    if not app.config['todo_include_todos']:
        for node in doctree.traverse(todo_node):
            node.parent.remove(node)

    # Replace all todolist nodes with a list of the collected todos.
    # Augment each todo with a backlink to the original location.
    env = app.builder.env

    if not hasattr(env, 'todo_all_todos'):
        env.todo_all_todos = []  # type: ignore

    for node in doctree.traverse(todolist):
        if not app.config['todo_include_todos']:
            node.replace_self([])
            continue

        content = []

        for todo_info in env.todo_all_todos:  # type: ignore
            para = nodes.paragraph(classes=['todo-source'])
            if app.config['todo_link_only']:
                description = _('<<original entry>>')
            else:
                description = (
                    _('(The <<original entry>> is located in %s, line %d.)') %
                    (todo_info['source'], todo_info['lineno'])
                )
            desc1 = description[:description.find('<<')]
            desc2 = description[description.find('>>') + 2:]
            para += nodes.Text(desc1, desc1)

            # Create a reference
            newnode = nodes.reference('', '', internal=True)
            innernode = nodes.emphasis(_('original entry'), _('original entry'))
            try:
                newnode['refuri'] = app.builder.get_relative_uri(fromdocname, todo_info['docname'])
                newnode['refuri'] += '#' + todo_info['target']['refid']
            # also trap KeyErrors
            #except NoUri
            except (NoUri, KeyError):
                # ignore if no URI can be determined, e.g. for LaTeX output
                pass
            newnode.append(innernode)
            para += newnode
            para += nodes.Text(desc2, desc2)

            # (Recursively) resolve references in the todo content
            todo_entry = todo_info['todo']
            env.resolve_references(todo_entry, todo_info['docname'], app.builder)

            # Insert into the todolist
            content.append(todo_entry)
            content.append(para)

        node.replace_self(content)
