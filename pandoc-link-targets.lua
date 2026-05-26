function Link(el)
  local target = el.target or ""

  if target:match("^#") or target:match("^mailto:") then
    return nil
  end

  el.attributes.target = "_blank"
  el.attributes.rel = "noopener"
  return el
end
