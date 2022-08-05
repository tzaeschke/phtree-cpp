// Copyright (c) Improbable Worlds Ltd, All Rights Reserved
#include "benchmark/logging.h"

namespace improbable::phtree::phbenchmark::logging {

void SetupDefaultLogging() {
    SetupLogging({}, spdlog::level::warn);
}

void SetupLogging(std::vector<spdlog::sink_ptr> sinks, spdlog::level::level_enum log_level) {
    auto& console_sink = sinks.emplace_back(std::make_shared<ConsoleSpdlogSink>());
    console_sink->set_level(log_level);

    // Find the minimum log level, in case one of the sinks passed to us has a lower log level.
    const auto& sink_with_lowest_log_level = *std::min_element(
        sinks.begin(),
        sinks.end(),
        [](const spdlog::sink_ptr& a, const spdlog::sink_ptr& b) -> bool {
            return a->level() < b->level();
        });
    spdlog::level::level_enum min_log_level =
        std::min(sink_with_lowest_log_level->level(), log_level);

    // Create the external logger, worker logger and the internal (default) logger from the same log
    // sinks. Each logsink can use `GetLoggerTypeFromMessage` to determine which logger a message
    // was logged to.
    spdlog::set_default_logger(
        std::make_shared<spdlog::logger>(kInternalLoggerName, sinks.begin(), sinks.end()));
    spdlog::set_level(min_log_level);
    spdlog::flush_on(min_log_level);
}

}  // namespace improbable::phtree::phbenchmark::logging
