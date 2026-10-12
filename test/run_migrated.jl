# Run the tests after migrating them to the new interface with tools/migration/migrate.jl, as a check of
# the migration tool: the migrated code (which is partially migrated where the tool cannot migrate it
# automatically) has to give the same results.
include(joinpath(@__DIR__, "..", "tools", "migration", "migrate.jl"))

function run_migrated_tests()
    srcfiles, messages = QAMigration.migrate_program(joinpath(@__DIR__, "runtests.jl"))
    println("Migrated the tests:")
    foreach(m -> println("  ", m), messages)
    dir = mktempdir()
    for sf in srcfiles
        out = joinpath(dir, relpath(sf.path, @__DIR__))
        mkpath(dirname(out))
        write(out, QAMigration.changed(sf) ? QAMigration.render(sf) : read(sf.path, String))
    end
    empty!(ARGS)
    include(joinpath(dir, "runtests.jl"))
end
run_migrated_tests()
